"""ChEMBL reference actives, Morgan novelty tiers, and PubChem PUG helpers."""

from __future__ import annotations

import json
import time
from collections import defaultdict
from pathlib import Path
from statistics import median
from typing import Any, Iterable

import pandas as pd
import requests
from rdkit import Chem

from hitpilot.fingerprints import fingerprint, max_tanimoto, morgan_generator
from hitpilot.io_utils import build_audit, load_compound_frame, write_audit, write_csv
from hitpilot.models import (
    ChemblActive,
    NoveltyConfig,
    NoveltyRecord,
    NoveltySummary,
    PubChemLookupResult,
    PubChemSimilarResult,
)
from hitpilot.version import __version__

CHEMBL_ACTIVITY_URL = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
PUBCHEM_PUG_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
DEFAULT_HEADERS = {
    "User-Agent": f"HitPilot/{__version__} (research; cheminformatics)",
    "Accept": "application/json",
}
NOVELTY_TIERS = (
    "known_or_close_match",
    "close_analog",
    "moderate_novelty",
    "novel_chemotype",
)


class ChemblFetchError(RuntimeError):
    """Raised when the ChEMBL REST API cannot be reached or parsed."""


def canonicalize_smiles(smiles: str) -> str | None:
    mol = Chem.MolFromSmiles("" if smiles is None else str(smiles).strip())
    if mol is None or mol.GetNumAtoms() == 0:
        return None
    return Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)


def novelty_tier(similarity: float, config: NoveltyConfig | None = None) -> str:
    config = config or NoveltyConfig()
    if similarity >= config.known_match_min:
        return "known_or_close_match"
    if similarity >= config.close_analog_min:
        return "close_analog"
    if similarity >= config.moderate_novelty_min:
        return "moderate_novelty"
    return "novel_chemotype"


def _cache_paths(cache_dir: str | Path, target_chembl_id: str) -> tuple[Path, Path]:
    cache_dir = Path(cache_dir)
    slug = target_chembl_id.strip().upper()
    return cache_dir / f"{slug}_actives.csv", cache_dir / f"{slug}_actives.meta.json"


def _write_chembl_cache(
    cache_csv: Path,
    cache_meta: Path,
    actives: list[ChemblActive],
    config: NoveltyConfig,
    n_activity_records: int,
) -> None:
    cache_csv.parent.mkdir(parents=True, exist_ok=True)
    rows = [
        {
            "molecule_chembl_id": item.molecule_chembl_id,
            "smiles": item.smiles,
            "canonical_smiles": item.canonical_smiles,
            "pchembl_value": item.pchembl_value,
            "target_chembl_id": item.target_chembl_id,
        }
        for item in actives
    ]
    pd.DataFrame(rows).to_csv(cache_csv, index=False)
    cache_meta.write_text(
        json.dumps(
            {
                "target_chembl_id": config.target_chembl_id,
                "pchembl_min": config.pchembl_min,
                "max_records": config.max_records,
                "n_activity_records": n_activity_records,
                "n_unique_actives": len(actives),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


def _load_chembl_cache(cache_csv: Path) -> list[ChemblActive]:
    frame = pd.read_csv(cache_csv)
    actives: list[ChemblActive] = []
    for row in frame.to_dict(orient="records"):
        smiles = str(row.get("canonical_smiles") or row.get("smiles") or "").strip()
        if not smiles:
            continue
        actives.append(
            ChemblActive(
                molecule_chembl_id=str(row["molecule_chembl_id"]),
                smiles=str(row.get("smiles") or smiles),
                canonical_smiles=smiles,
                pchembl_value=float(row["pchembl_value"]),
                target_chembl_id=str(row.get("target_chembl_id") or ""),
            )
        )
    return actives


def _cache_matches(cache_meta: Path, config: NoveltyConfig) -> bool:
    if not cache_meta.is_file():
        return False
    try:
        meta = json.loads(cache_meta.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return False
    return (
        str(meta.get("target_chembl_id", "")).upper() == config.target_chembl_id.upper()
        and float(meta.get("pchembl_min", -1)) == float(config.pchembl_min)
        and int(meta.get("max_records", -1)) == int(config.max_records)
    )


def _deduplicate_actives(
    rows: Iterable[dict[str, Any]],
    target_chembl_id: str,
) -> list[ChemblActive]:
    best: dict[str, ChemblActive] = {}
    for row in rows:
        raw = row.get("canonical_smiles") or row.get("smiles")
        canon = canonicalize_smiles(str(raw) if raw is not None else "")
        if canon is None:
            continue
        mol_id = str(row.get("molecule_chembl_id") or "").strip()
        try:
            pchembl = float(row.get("pchembl_value"))
        except (TypeError, ValueError):
            continue
        current = best.get(canon)
        if current is None or pchembl > current.pchembl_value:
            best[canon] = ChemblActive(
                molecule_chembl_id=mol_id or "UNKNOWN",
                smiles=str(raw).strip(),
                canonical_smiles=canon,
                pchembl_value=pchembl,
                target_chembl_id=target_chembl_id,
            )
    return sorted(best.values(), key=lambda item: (-item.pchembl_value, item.molecule_chembl_id))


def fetch_chembl_activities(
    *,
    target_chembl_id: str,
    pchembl_min: float,
    max_records: int,
    timeout: int = 15,
    session: requests.Session | None = None,
    page_size: int = 1000,
) -> list[dict[str, Any]]:
    """Page the live ChEMBL activity API, most potent first, up to *max_records*."""
    client = session or requests.Session()
    collected: list[dict[str, Any]] = []
    offset = 0
    page_size = min(max(1, page_size), max(1, max_records), 1000)
    while len(collected) < max_records:
        limit = min(page_size, max_records - len(collected))
        params = {
            "target_chembl_id": target_chembl_id,
            "pchembl_value__gte": pchembl_min,
            "limit": limit,
            "offset": offset,
            "order_by": "-pchembl_value",
        }
        try:
            response = client.get(
                CHEMBL_ACTIVITY_URL,
                params=params,
                timeout=(timeout, max(timeout, 60)),
                headers=DEFAULT_HEADERS,
            )
            response.raise_for_status()
            payload = response.json()
        except Exception as exc:  # noqa: BLE001
            raise ChemblFetchError(f"ChEMBL activity request failed: {exc}") from exc
        activities = payload.get("activities") or []
        if not activities:
            break
        collected.extend(activities)
        offset += len(activities)
        total = (payload.get("page_meta") or {}).get("total_count")
        if total is not None and offset >= int(total):
            break
        if len(activities) < limit:
            break
        time.sleep(0.15)
    return collected[:max_records]


def fetch_chembl_target_actives(
    target_chembl_id: str = "CHEMBL206",
    pchembl_min: float = 6.0,
    max_records: int = 2000,
    cache_dir: str | Path = "data/chembl",
    *,
    refresh_cache: bool = False,
    timeout: int = 15,
    session: requests.Session | None = None,
    page_size: int = 1000,
) -> list[ChemblActive]:
    """Fetch unique human ChEMBL actives and cache them under data/chembl/."""
    config = NoveltyConfig(
        target_chembl_id=target_chembl_id,
        pchembl_min=pchembl_min,
        max_records=max_records,
        cache_dir=str(cache_dir),
        refresh_cache=refresh_cache,
        timeout=timeout,
    )
    cache_csv, cache_meta = _cache_paths(cache_dir, target_chembl_id)
    if not refresh_cache and cache_csv.is_file() and _cache_matches(cache_meta, config):
        return _load_chembl_cache(cache_csv)

    raw = fetch_chembl_activities(
        target_chembl_id=target_chembl_id,
        pchembl_min=pchembl_min,
        max_records=max_records,
        timeout=timeout,
        session=session,
        page_size=page_size,
    )
    actives = _deduplicate_actives(raw, target_chembl_id)
    _write_chembl_cache(cache_csv, cache_meta, actives, config, n_activity_records=len(raw))
    return actives


def _parse_pubchem_cids(payload: Any) -> list[int]:
    if not isinstance(payload, dict):
        return []
    identifier = payload.get("IdentifierList") or payload.get("identifierList") or {}
    cids = identifier.get("CID") or identifier.get("cid") or []
    if isinstance(cids, int):
        return [cids]
    result: list[int] = []
    for value in cids:
        try:
            result.append(int(value))
        except (TypeError, ValueError):
            continue
    return result


def _pubchem_sleep(min_interval: float, last_call: list[float]) -> None:
    if last_call:
        elapsed = time.monotonic() - last_call[0]
        if elapsed < min_interval:
            time.sleep(min_interval - elapsed)
    last_call[:] = [time.monotonic()]


_PUBCHEM_LAST_CALL: list[float] = []


def lookup_pubchem_compound(
    smiles: str,
    *,
    timeout: int = 15,
    min_interval: float = 0.25,
    session: requests.Session | None = None,
) -> PubChemLookupResult:
    """Look up PubChem CID(s) for one SMILES (rate-limited PUG REST)."""
    query = "" if smiles is None else str(smiles).strip()
    if not query:
        return PubChemLookupResult(query_smiles="", error="empty SMILES")
    _pubchem_sleep(min_interval, _PUBCHEM_LAST_CALL)
    client = session or requests
    url = f"{PUBCHEM_PUG_BASE}/compound/smiles/cids/JSON"
    try:
        response = client.post(
            url,
            data={"smiles": query},
            timeout=(timeout, max(timeout, 30)),
            headers=DEFAULT_HEADERS,
        )
        if response.status_code == 404:
            return PubChemLookupResult(query_smiles=query, cids=[])
        response.raise_for_status()
        cids = _parse_pubchem_cids(response.json())
        return PubChemLookupResult(query_smiles=query, cids=cids)
    except Exception as exc:  # noqa: BLE001
        return PubChemLookupResult(query_smiles=query, error=str(exc))


def search_pubchem_similar(
    smiles: str,
    threshold: int = 85,
    max_records: int = 5,
    *,
    timeout: int = 15,
    min_interval: float = 0.25,
    session: requests.Session | None = None,
) -> PubChemSimilarResult:
    """2D similarity search for a shortlist candidate (rate-limited PUG REST)."""
    query = "" if smiles is None else str(smiles).strip()
    if not query:
        return PubChemSimilarResult(
            query_smiles="",
            threshold=threshold,
            max_records=max_records,
            error="empty SMILES",
        )
    _pubchem_sleep(min_interval, _PUBCHEM_LAST_CALL)
    client = session or requests
    url = f"{PUBCHEM_PUG_BASE}/compound/fastsimilarity_2d/smiles/cids/JSON"
    try:
        response = client.post(
            url,
            data={"smiles": query},
            params={"Threshold": int(threshold), "MaxRecords": int(max_records)},
            timeout=(timeout, max(timeout, 30)),
            headers=DEFAULT_HEADERS,
        )
        if response.status_code == 404:
            return PubChemSimilarResult(
                query_smiles=query,
                threshold=threshold,
                max_records=max_records,
                cids=[],
            )
        response.raise_for_status()
        cids = _parse_pubchem_cids(response.json())[: int(max_records)]
        return PubChemSimilarResult(
            query_smiles=query,
            threshold=threshold,
            max_records=max_records,
            cids=cids,
        )
    except Exception as exc:  # noqa: BLE001
        return PubChemSimilarResult(
            query_smiles=query,
            threshold=threshold,
            max_records=max_records,
            error=str(exc),
        )


def run_novelty(
    rows: Iterable[dict[str, Any]],
    references: list[ChemblActive],
    config: NoveltyConfig | None = None,
    *,
    smiles_field: str = "smiles",
    id_field: str = "compound_id",
    label_field: str = "label",
) -> tuple[list[NoveltyRecord], NoveltySummary]:
    config = config or NoveltyConfig()
    generator = morgan_generator(radius=config.radius, fp_size=config.fp_size)
    ref_mols: list[ChemblActive] = []
    ref_fps = []
    for item in references:
        mol = Chem.MolFromSmiles(item.canonical_smiles)
        if mol is None:
            continue
        ref_mols.append(item)
        ref_fps.append(fingerprint(mol, generator))

    records: list[NoveltyRecord] = []
    for index, row in enumerate(rows):
        compound_id = str(row.get(id_field) or row.get("compound_id") or f"cmpd_{index + 1:06d}")
        smiles = str(row.get(smiles_field, row.get("smiles", "")) or "")
        raw_label = row.get(label_field, row.get("label"))
        label = None if raw_label is None or str(raw_label).strip() == "" else str(raw_label)
        mol = Chem.MolFromSmiles(smiles.strip()) if smiles.strip() else None
        if mol is None:
            records.append(
                NoveltyRecord(
                    compound_id=compound_id,
                    smiles=smiles,
                    label=label,
                    error="invalid SMILES" if smiles.strip() else "empty SMILES",
                )
            )
            continue
        fp = fingerprint(mol, generator)
        sim, idx = max_tanimoto(fp, ref_fps)
        if idx < 0:
            records.append(
                NoveltyRecord(
                    compound_id=compound_id,
                    smiles=smiles,
                    max_tanimoto=0.0,
                    novelty_tier="novel_chemotype",
                    label=label,
                )
            )
            continue
        nearest = ref_mols[idx]
        records.append(
            NoveltyRecord(
                compound_id=compound_id,
                smiles=smiles,
                max_tanimoto=sim,
                nearest_ref_id=nearest.molecule_chembl_id,
                nearest_ref_pchembl=nearest.pchembl_value,
                nearest_ref_smiles=nearest.canonical_smiles,
                novelty_tier=novelty_tier(sim, config),
                label=label,
            )
        )
    return records, summarize_novelty(records, config, n_reference_actives=len(ref_mols))


def summarize_novelty(
    records: list[NoveltyRecord],
    config: NoveltyConfig,
    *,
    n_reference_actives: int,
) -> NoveltySummary:
    valid = [rec for rec in records if rec.error is None and rec.max_tanimoto is not None]
    counts: dict[str, int] = {tier: 0 for tier in NOVELTY_TIERS}
    for rec in valid:
        if rec.novelty_tier:
            counts[rec.novelty_tier] = counts.get(rec.novelty_tier, 0) + 1
    by_label: dict[str, dict[str, int]] = {}
    grouped: dict[str, list[NoveltyRecord]] = defaultdict(list)
    for rec in records:
        if rec.label:
            grouped[str(rec.label)].append(rec)
    for label, subset in sorted(grouped.items()):
        sub_counts = {tier: 0 for tier in NOVELTY_TIERS}
        sub_valid = [rec for rec in subset if rec.error is None]
        for rec in sub_valid:
            if rec.novelty_tier:
                sub_counts[rec.novelty_tier] = sub_counts.get(rec.novelty_tier, 0) + 1
        by_label[label] = {"n_total": len(subset), "n_valid": len(sub_valid), **sub_counts}
    sims = [float(rec.max_tanimoto) for rec in valid if rec.max_tanimoto is not None]
    return NoveltySummary(
        n_total=len(records),
        n_valid=len(valid),
        n_invalid=len(records) - len(valid),
        n_reference_actives=n_reference_actives,
        target_chembl_id=config.target_chembl_id,
        pchembl_min=config.pchembl_min,
        median_max_tanimoto=float(median(sims)) if sims else None,
        counts_by_tier=counts,
        by_label=by_label,
    )


def run_novelty_file(
    input_path: str | Path,
    output_dir: str | Path,
    config: NoveltyConfig | None = None,
    *,
    smiles_column: str | None = None,
    id_column: str | None = None,
    label_column: str | None = None,
    references: list[ChemblActive] | None = None,
):
    """Score novelty against ChEMBL reference actives and write CSV plus audit."""
    config = config or NoveltyConfig()
    input_path = Path(input_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    if references is None:
        references = fetch_chembl_target_actives(
            target_chembl_id=config.target_chembl_id,
            pchembl_min=config.pchembl_min,
            max_records=config.max_records,
            cache_dir=config.cache_dir,
            refresh_cache=config.refresh_cache,
            timeout=config.timeout,
        )
    frame = load_compound_frame(
        input_path,
        smiles_column=smiles_column,
        id_column=id_column,
        label_column=label_column,
    )
    records, summary = run_novelty(frame.to_dict(orient="records"), references, config)
    csv_path = output_dir / "novelty.csv"
    write_csv(csv_path, [rec.to_csv_row() for rec in records])
    manifest = build_audit(
        command="novelty",
        config=config,
        input_path=input_path,
        output_dir=output_dir,
        summary=summary,
        outputs=[str(csv_path.resolve()), str((output_dir / "audit_log.json").resolve())],
    )
    write_audit(output_dir / "audit_log.json", manifest)
    return records, summary, manifest
