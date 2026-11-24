"""Unit tests for RDKit hygiene, alerts, and audit logs."""

from __future__ import annotations

import json
from pathlib import Path

from hitpilot.hygiene import (
    DROP_BRENK,
    DROP_DUPLICATE,
    DROP_EMPTY,
    DROP_INORGANIC,
    DROP_INVALID,
    DROP_NIH,
    DROP_PAINS,
    DROP_SANITIZE,
    run_hygiene,
    run_hygiene_file,
    standardize_smiles,
)
from hitpilot.models import HygieneConfig

# Public, well-known structures (not invented identifiers).
ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
IBUPROFEN = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
ERLOTINIB = "C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1"
RHODANINE_BENZYLIDENE = "O=C1NC(=S)SC1=Cc1ccccc1"
CATECHOL = "Oc1ccccc1O"
QUINONE = "O=C1C=CC(=O)C=C1"
ACYL_CHLORIDE = "CC(=O)Cl"
EPOXIDE = "C1CO1"
NITROBENZENE = "O=[N+]([O-])c1ccccc1"


def _by_id(records):
    return {rec.compound_id: rec for rec in records}


def test_invalid_smiles_is_captured():
    result = standardize_smiles("not_a_smiles")
    assert result.mol is None
    assert result.drop_reasons == [DROP_INVALID]
    assert result.error


def test_empty_smiles_is_invalid():
    result = standardize_smiles("   ")
    assert result.drop_reasons == [DROP_INVALID]


def test_impossible_valence_fails_sanitization():
    # Pentavalent carbon is a real RDKit AtomValenceException.
    result = standardize_smiles("C(C)(C)(C)(C)C")
    assert result.drop_reasons == [DROP_SANITIZE]
    assert result.error
    records, summary = run_hygiene(
        [{"compound_id": "bad", "smiles": "C(C)(C)(C)(C)C"}],
        HygieneConfig(),
    )
    assert records[0].drop_reasons == [DROP_SANITIZE]
    assert records[0].passed is False
    assert summary.counts_by_drop_reason[DROP_SANITIZE] == 1


def test_salt_counterion_and_solvent_stripping():
    cases = {
        "amine_hcl": "CCNCC.Cl",
        "tfa_salt": "CCNCC.OC(=O)C(F)(F)F",
        "with_water": "CCNCC.O",
        "free_base": "CCNCC",
    }
    keys = {}
    for name, smiles in cases.items():
        result = standardize_smiles(smiles)
        assert result.drop_reasons == [], result.error
        assert result.smiles is not None
        keys[name] = (result.smiles, result.inchikey)
        assert "Cl" not in result.smiles
        assert "F" not in result.smiles
    assert keys["amine_hcl"] == keys["free_base"] == keys["tfa_salt"] == keys["with_water"]
    assert keys["free_base"][1] == "HPNMFZURTQLUMO-UHFFFAOYSA-N"


def test_sodium_salt_is_neutralized():
    result = standardize_smiles("CC(=O)[O-].[Na+]")
    assert result.drop_reasons == []
    assert result.smiles == "CC(=O)O"
    assert result.inchikey == "QTBSBXVTEAMEQO-UHFFFAOYSA-N"


def test_canonicalization_of_equivalent_smiles():
    aromatic = standardize_smiles("c1ccccc1")
    kekule = standardize_smiles("C1=CC=CC=C1")
    assert aromatic.smiles == kekule.smiles == "c1ccccc1"
    assert aromatic.inchikey == kekule.inchikey


def test_inorganic_and_empty_structures_are_rejected():
    water = standardize_smiles("O")
    assert DROP_INORGANIC in water.drop_reasons
    salt = standardize_smiles("[Na+].[Cl-]")
    assert DROP_INORGANIC in salt.drop_reasons or DROP_EMPTY in salt.drop_reasons
    sodium = standardize_smiles("[Na+]")
    assert DROP_INORGANIC in sodium.drop_reasons


def test_deduplication_of_equivalent_and_salted_structures():
    rows = [
        {"compound_id": "benzene_arom", "smiles": "c1ccccc1"},
        {"compound_id": "benzene_kekule", "smiles": "C1=CC=CC=C1"},
        {"compound_id": "diethylamine", "smiles": "CCNCC"},
        {"compound_id": "diethylamine_hcl", "smiles": "CCNCC.Cl"},
    ]
    records, summary = run_hygiene(rows, HygieneConfig())
    recs = _by_id(records)
    assert recs["benzene_arom"].passed is True
    assert recs["benzene_kekule"].drop_reasons == [DROP_DUPLICATE]
    assert recs["benzene_kekule"].duplicate_of == "benzene_arom"
    assert recs["diethylamine"].passed is True
    assert recs["diethylamine_hcl"].drop_reasons == [DROP_DUPLICATE]
    assert recs["diethylamine_hcl"].duplicate_of == "diethylamine"
    assert summary.counts_by_drop_reason[DROP_DUPLICATE] == 2
    assert summary.n_passed == 2


def test_score_aware_deduplication_keeps_best():
    rows = [
        {"compound_id": "low", "smiles": "CCO", "score": 1.2},
        {"compound_id": "high", "smiles": "OCC", "score": 5.0},
    ]
    descending, _ = run_hygiene(rows, HygieneConfig(score_column="score", score_ascending=False))
    recs = _by_id(descending)
    assert recs["high"].passed is True
    assert recs["low"].duplicate_of == "high"
    assert DROP_DUPLICATE in recs["low"].drop_reasons

    ascending, _ = run_hygiene(rows, HygieneConfig(score_column="score", score_ascending=True))
    recs = _by_id(ascending)
    assert recs["low"].passed is True
    assert recs["high"].duplicate_of == "low"


def test_real_pains_triggers_rhodanine_catechol_quinone():
    rows = [
        {"compound_id": "rhodanine", "smiles": RHODANINE_BENZYLIDENE},
        {"compound_id": "catechol", "smiles": CATECHOL},
        {"compound_id": "quinone", "smiles": QUINONE},
    ]
    records, _ = run_hygiene(rows, HygieneConfig())
    recs = _by_id(records)

    rhod_catalogs = {hit.catalog for hit in recs["rhodanine"].alerts}
    rhod_rules = {hit.rule for hit in recs["rhodanine"].alerts}
    assert "PAINS_A" in rhod_catalogs
    assert "ene_rhod_A(235)" in rhod_rules
    assert "NIH" in rhod_catalogs and "rhodanine" in rhod_rules
    assert DROP_PAINS in recs["rhodanine"].drop_reasons
    assert recs["rhodanine"].passed is False

    cat_rules = {hit.rule for hit in recs["catechol"].alerts}
    assert "catechol_A(92)" in cat_rules
    assert DROP_PAINS in recs["catechol"].drop_reasons

    quin_rules = {hit.rule for hit in recs["quinone"].alerts}
    quin_catalogs = {hit.catalog for hit in recs["quinone"].alerts}
    assert "quinone_A(370)" in quin_rules
    assert "PAINS_A" in quin_catalogs
    assert "para_quinone" in quin_rules or "bis_keto_olefin" in quin_rules


def test_real_brenk_triggers_acyl_halide_epoxide_nitro():
    rows = [
        {"compound_id": "acyl", "smiles": ACYL_CHLORIDE},
        {"compound_id": "epoxide", "smiles": EPOXIDE},
        {"compound_id": "nitro", "smiles": NITROBENZENE},
    ]
    records, _ = run_hygiene(rows, HygieneConfig())
    recs = _by_id(records)

    assert any(hit.catalog == "BRENK" and hit.rule == "acid_halide" for hit in recs["acyl"].alerts)
    assert DROP_BRENK in recs["acyl"].drop_reasons
    # NIH also annotates acyl halides as carbonyl_halide.
    assert any(
        hit.catalog == "NIH" and hit.rule == "carbonyl_halide"
        for hit in recs["acyl"].alerts
    )
    assert DROP_NIH in recs["acyl"].drop_reasons

    assert any(
        hit.catalog == "BRENK" and hit.rule == "Three-membered_heterocycle"
        for hit in recs["epoxide"].alerts
    )
    assert DROP_BRENK in recs["epoxide"].drop_reasons

    nitro_rules = {hit.rule for hit in recs["nitro"].alerts if hit.catalog == "BRENK"}
    assert "nitro_group" in nitro_rules
    assert DROP_BRENK in recs["nitro"].drop_reasons
    assert recs["nitro"].passed is False


def test_nih_triggers_and_flag_only_mode():
    records, _ = run_hygiene(
        [{"compound_id": "aldehyde", "smiles": "O=Cc1ccccc1"}],
        HygieneConfig(drop_brenk=False, drop_nih=False),
    )
    rec = records[0]
    catalogs = {hit.catalog for hit in rec.alerts}
    rules = {hit.rule for hit in rec.alerts}
    assert "NIH" in catalogs
    assert "aldehyde" in rules
    assert rec.passed is True
    assert DROP_NIH not in rec.drop_reasons
    assert DROP_BRENK not in rec.drop_reasons


def test_known_drugs_against_real_filtercatalog_behavior():
    rows = [
        {"compound_id": "aspirin", "smiles": ASPIRIN},
        {"compound_id": "ibuprofen", "smiles": IBUPROFEN},
        {"compound_id": "erlotinib", "smiles": ERLOTINIB},
    ]
    records, _ = run_hygiene(rows, HygieneConfig())
    recs = _by_id(records)

    # Ibuprofen has no PAINS, Brenk, or NIH alerts in the RDKit catalogs.
    assert recs["ibuprofen"].alerts == []
    assert recs["ibuprofen"].passed is True
    assert recs["ibuprofen"].inchikey == "HEFNNWSXXWATRW-UHFFFAOYSA-N"

    # Aspirin is not a PAINS/NIH hit; Brenk flags the phenol ester.
    aspirin_catalogs = {hit.catalog for hit in recs["aspirin"].alerts}
    assert "PAINS_A" not in aspirin_catalogs
    assert "PAINS_B" not in aspirin_catalogs
    assert "PAINS_C" not in aspirin_catalogs
    assert "NIH" not in aspirin_catalogs
    assert any(
        hit.catalog == "BRENK" and hit.rule == "phenol_ester"
        for hit in recs["aspirin"].alerts
    )
    assert recs["aspirin"].passed is False

    # Erlotinib is not a PAINS/NIH hit; Brenk flags the alkyne / long chain.
    erlotinib_catalogs = {hit.catalog for hit in recs["erlotinib"].alerts}
    assert not any(name.startswith("PAINS") for name in erlotinib_catalogs)
    assert "NIH" not in erlotinib_catalogs
    assert "BRENK" in erlotinib_catalogs
    assert recs["erlotinib"].passed is False

    # Flag-only Brenk keeps both drugs that would otherwise be dropped.
    kept, _ = run_hygiene(rows, HygieneConfig(drop_brenk=False))
    kept_recs = _by_id(kept)
    assert kept_recs["aspirin"].passed is True
    assert kept_recs["erlotinib"].passed is True
    assert kept_recs["ibuprofen"].passed is True


def test_audit_log_json_integrity(tmp_path: Path):
    input_path = tmp_path / "input.csv"
    input_path.write_text(
        "compound_id,smiles,label,score\n"
        f"ibu,{IBUPROFEN},active,9.1\n"
        f"cat,{CATECHOL},inactive,0.2\n"
        "dup,CCO,active,1.0\n"
        "dup2,OCC,active,0.4\n",
        encoding="utf-8",
    )
    output_dir = tmp_path / "out"
    manifest = run_hygiene_file(
        input_path,
        output_dir,
        HygieneConfig(score_column="score", score_ascending=False),
        score_column="score",
    )
    audit_path = output_dir / "audit_log.json"
    all_path = output_dir / "all_compounds.csv"
    passed_path = output_dir / "passed_compounds.csv"
    assert audit_path.is_file() and audit_path.stat().st_size > 0
    assert all_path.is_file() and all_path.stat().st_size > 0
    assert passed_path.is_file() and passed_path.stat().st_size > 0

    payload = json.loads(audit_path.read_text(encoding="utf-8"))
    assert payload["hitpilot_version"] == "0.1.0"
    assert payload["rdkit_version"]
    assert payload["timestamp_utc"].endswith("Z")
    assert payload["input_sha256"] == manifest.input_sha256
    assert len(payload["input_sha256"]) == 64
    assert payload["summary"]["n_total"] == 4
    assert "pains_alert" in payload["summary"]["counts_by_drop_reason"]
    assert "duplicate_structure" in payload["summary"]["counts_by_drop_reason"]
    assert "PAINS_B" in payload["summary"]["counts_by_catalog"]
    assert payload["config"]["drop_pains"] is True
    assert payload["summary"]["by_label"]["active"]["n_total"] == 3
