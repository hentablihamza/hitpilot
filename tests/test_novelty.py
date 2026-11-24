"""Unit tests for ChEMBL fetching, novelty tiers, and PubChem helpers."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from hitpilot.models import ChemblActive, NoveltyConfig
from hitpilot.novelty import (
    _deduplicate_actives,
    canonicalize_smiles,
    fetch_chembl_target_actives,
    lookup_pubchem_compound,
    novelty_tier,
    run_novelty,
    run_novelty_file,
    search_pubchem_similar,
)

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
IBUPROFEN = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
ETHANOL = "CCO"


class FakeResponse:
    def __init__(self, payload, status_code=200):
        self._payload = payload
        self.status_code = status_code
        self.content = json.dumps(payload).encode()

    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(f"HTTP {self.status_code}")

    def json(self):
        return self._payload


class FakeSession:
    def __init__(self, handler):
        self.handler = handler
        self.calls: list[dict] = []

    def get(self, url, params=None, timeout=None, headers=None):
        self.calls.append({"url": url, "params": params})
        return self.handler("GET", url, params=params)

    def post(self, url, data=None, params=None, timeout=None, headers=None):
        self.calls.append({"url": url, "data": data, "params": params})
        return self.handler("POST", url, data=data, params=params)


def test_novelty_tier_boundaries():
    config = NoveltyConfig()
    assert novelty_tier(0.85, config) == "known_or_close_match"
    assert novelty_tier(0.99, config) == "known_or_close_match"
    assert novelty_tier(0.65, config) == "close_analog"
    assert novelty_tier(0.84, config) == "close_analog"
    assert novelty_tier(0.45, config) == "moderate_novelty"
    assert novelty_tier(0.64, config) == "moderate_novelty"
    assert novelty_tier(0.449, config) == "novel_chemotype"
    assert novelty_tier(0.0, config) == "novel_chemotype"


def test_dedup_keeps_highest_pchembl_and_canonical_smiles():
    rows = [
        {"molecule_chembl_id": "CHEMBL_LOW", "canonical_smiles": "CCO", "pchembl_value": "6.2"},
        {"molecule_chembl_id": "CHEMBL_HIGH", "canonical_smiles": "OCC", "pchembl_value": "9.1"},
        {"molecule_chembl_id": "CHEMBL25", "canonical_smiles": ASPIRIN, "pchembl_value": "7.4"},
        {"molecule_chembl_id": "BAD", "canonical_smiles": "not_a_smiles", "pchembl_value": "9.9"},
    ]
    actives = _deduplicate_actives(rows, "CHEMBL206")
    by_id = {item.molecule_chembl_id: item for item in actives}
    assert set(by_id) == {"CHEMBL_HIGH", "CHEMBL25"}
    assert by_id["CHEMBL_HIGH"].canonical_smiles == canonicalize_smiles("CCO")
    assert by_id["CHEMBL_HIGH"].pchembl_value == 9.1


def test_fetch_chembl_paginates_and_caches(tmp_path: Path):
    pages = {
        0: {
            "activities": [
                {
                    "molecule_chembl_id": "CHEMBL521",
                    "canonical_smiles": IBUPROFEN,
                    "pchembl_value": "7.20",
                },
                {
                    "molecule_chembl_id": "CHEMBL25",
                    "canonical_smiles": ASPIRIN,
                    "pchembl_value": "6.10",
                },
            ],
            "page_meta": {"total_count": 3, "limit": 2, "offset": 0},
        },
        2: {
            "activities": [
                {
                    "molecule_chembl_id": "CHEMBL521B",
                    "canonical_smiles": IBUPROFEN,  # same structure, lower pChEMBL
                    "pchembl_value": "6.05",
                }
            ],
            "page_meta": {"total_count": 3, "limit": 2, "offset": 2},
        },
    }

    def handler(method, url, params=None, **kwargs):
        assert "activity.json" in url
        return FakeResponse(pages[params["offset"]])

    session = FakeSession(handler)
    actives = fetch_chembl_target_actives(
        "CHEMBL206",
        pchembl_min=6.0,
        max_records=3,
        cache_dir=tmp_path,
        refresh_cache=True,
        session=session,
        page_size=2,
    )
    assert len(session.calls) == 2
    assert {item.molecule_chembl_id for item in actives} == {"CHEMBL521", "CHEMBL25"}
    ibu = next(item for item in actives if item.molecule_chembl_id == "CHEMBL521")
    assert ibu.pchembl_value == 7.20
    cache = tmp_path / "CHEMBL206_actives.csv"
    assert cache.is_file() and cache.stat().st_size > 0

    # Second call should hit the cache and not touch the session.
    session.calls.clear()
    cached = fetch_chembl_target_actives(
        "CHEMBL206",
        pchembl_min=6.0,
        max_records=3,
        cache_dir=tmp_path,
        refresh_cache=False,
        session=session,
    )
    assert session.calls == []
    assert {item.molecule_chembl_id for item in cached} == {"CHEMBL521", "CHEMBL25"}


def test_nearest_neighbor_novelty_on_known_drugs():
    refs = [
        ChemblActive(
            molecule_chembl_id="CHEMBL521",
            smiles=IBUPROFEN,
            canonical_smiles=canonicalize_smiles(IBUPROFEN),
            pchembl_value=7.2,
            target_chembl_id="CHEMBL206",
        ),
        ChemblActive(
            molecule_chembl_id="CHEMBL25",
            smiles=ASPIRIN,
            canonical_smiles=canonicalize_smiles(ASPIRIN),
            pchembl_value=6.1,
            target_chembl_id="CHEMBL206",
        ),
    ]
    rows = [
        {"compound_id": "self", "smiles": IBUPROFEN, "label": "active"},
        {"compound_id": "ethanol", "smiles": ETHANOL, "label": "inactive"},
        {"compound_id": "bad", "smiles": "not_a_smiles", "label": "inactive"},
    ]
    records, summary = run_novelty(rows, refs, NoveltyConfig())
    recs = {rec.compound_id: rec for rec in records}
    assert recs["self"].max_tanimoto == pytest.approx(1.0)
    assert recs["self"].nearest_ref_id == "CHEMBL521"
    assert recs["self"].nearest_ref_pchembl == 7.2
    assert recs["self"].novelty_tier == "known_or_close_match"
    assert recs["ethanol"].max_tanimoto < 0.45
    assert recs["ethanol"].novelty_tier == "novel_chemotype"
    assert recs["bad"].error == "invalid SMILES"
    assert summary.n_reference_actives == 2
    assert summary.counts_by_tier["known_or_close_match"] == 1
    assert summary.counts_by_tier["novel_chemotype"] == 1
    assert summary.by_label["active"]["known_or_close_match"] == 1


def test_pubchem_helpers_are_rate_limited_and_parse_cids():
    def handler(method, url, data=None, params=None, **kwargs):
        assert method == "POST"
        if "fastsimilarity_2d" in url:
            assert params["Threshold"] == 85
            assert params["MaxRecords"] == 5
            return FakeResponse({"IdentifierList": {"CID": [2244, 3672, 135]}})
        return FakeResponse({"IdentifierList": {"CID": [2244]}})

    session = FakeSession(handler)
    identity = lookup_pubchem_compound(ASPIRIN, min_interval=0.0, session=session)
    assert identity.cids == [2244]
    assert identity.error is None
    similar = search_pubchem_similar(
        ASPIRIN, threshold=85, max_records=5, min_interval=0.0, session=session
    )
    assert similar.cids == [2244, 3672, 135]
    assert similar.threshold == 85


def test_novelty_file_writes_audit(tmp_path: Path):
    input_path = tmp_path / "in.csv"
    input_path.write_text(f"compound_id,smiles\nibu,{IBUPROFEN}\n", encoding="utf-8")
    refs = [
        ChemblActive(
            molecule_chembl_id="CHEMBL521",
            smiles=IBUPROFEN,
            canonical_smiles=canonicalize_smiles(IBUPROFEN),
            pchembl_value=7.2,
            target_chembl_id="CHEMBL206",
        )
    ]
    records, summary, manifest = run_novelty_file(
        input_path, tmp_path / "out", NoveltyConfig(), references=refs
    )
    assert records[0].novelty_tier == "known_or_close_match"
    audit = json.loads((tmp_path / "out" / "audit_log.json").read_text(encoding="utf-8"))
    assert audit["command"] == "novelty"
    assert audit["summary"]["n_reference_actives"] == 1
    assert manifest.input_sha256 == audit["input_sha256"]
    assert summary.n_total == 1
