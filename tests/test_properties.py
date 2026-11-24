"""Unit tests for RDKit descriptors, Ertl SA score, and Lipinski/Veber flags."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from rdkit import Chem, RDConfig

from hitpilot.models import PropertiesConfig
from hitpilot.properties import (
    compute_properties,
    load_sascorer,
    property_record,
    run_properties,
    run_properties_file,
    sa_score,
)

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
IBUPROFEN = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
BENZENE = "c1ccccc1"
STEARIC = "CCCCCCCCCCCCCCCCCC(=O)O"
PEG12 = "CCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCO"


def test_sascorer_loads_from_rdkit_contrib():
    module = load_sascorer()
    expected = Path(RDConfig.RDContribDir) / "SA_Score" / "sascorer.py"
    assert expected.is_file()
    assert module.__file__ == str(expected)
    assert hasattr(module, "calculateScore")


def test_known_drug_descriptors_and_sa_scores():
    aspirin = compute_properties(Chem.MolFromSmiles(ASPIRIN))
    assert aspirin["mw"] == pytest.approx(180.15899999999996)
    assert aspirin["clogp"] == pytest.approx(1.3101)
    assert aspirin["tpsa"] == pytest.approx(63.60000000000001)
    assert aspirin["hbd"] == 1
    assert aspirin["hba"] == 4
    assert aspirin["rotatable_bonds"] == 2
    assert aspirin["fraction_csp3"] == pytest.approx(0.1111111111111111)
    assert aspirin["ring_count"] == 1
    assert aspirin["qed"] == pytest.approx(0.5501217966938848)
    assert aspirin["sa_score"] == pytest.approx(1.580039750008826)
    assert aspirin["lipinski_violations"] == 0
    assert aspirin["lipinski_pass"] is True
    assert aspirin["lipinski_strict_pass"] is True
    assert aspirin["veber_pass"] is True

    ibuprofen = compute_properties(Chem.MolFromSmiles(IBUPROFEN))
    assert ibuprofen["mw"] == pytest.approx(206.28499999999997)
    assert ibuprofen["clogp"] == pytest.approx(3.073200000000001)
    assert ibuprofen["hbd"] == 1 and ibuprofen["hba"] == 2
    assert ibuprofen["sa_score"] == pytest.approx(2.1917554204925054)
    assert ibuprofen["lipinski_strict_pass"] is True
    assert ibuprofen["veber_pass"] is True

    benzene = compute_properties(Chem.MolFromSmiles(BENZENE))
    assert benzene["sa_score"] == pytest.approx(1.0)
    assert benzene["ring_count"] == 1
    assert sa_score(Chem.MolFromSmiles(BENZENE)) == pytest.approx(benzene["sa_score"])


def test_lipinski_and_veber_thresholds():
    stearic = compute_properties(Chem.MolFromSmiles(STEARIC))
    # cLogP > 5, one Lipinski violation; 16 rotatable bonds fail Veber.
    assert stearic["clogp"] > 5
    assert stearic["mw"] < 500
    assert stearic["lipinski_violations"] == 1
    assert stearic["lipinski_pass"] is True
    assert stearic["lipinski_strict_pass"] is False
    assert stearic["rotatable_bonds"] > 10
    assert stearic["veber_pass"] is False

    peg = compute_properties(Chem.MolFromSmiles(PEG12))
    # MW > 500 and HBA > 10: two violations, Lipinski fail; many rotors fail Veber.
    assert peg["mw"] > 500
    assert peg["hba"] > 10
    assert peg["lipinski_violations"] == 2
    assert peg["lipinski_pass"] is False
    assert peg["lipinski_strict_pass"] is False
    assert peg["veber_pass"] is False


def test_invalid_smiles_are_captured():
    rec = property_record("bad", "not_a_smiles")
    assert rec.error == "invalid SMILES"
    assert rec.mw is None
    assert rec.lipinski_pass is None
    records, summary = run_properties(
        [
            {"compound_id": "bad", "smiles": "not_a_smiles"},
            {"compound_id": "ibu", "smiles": IBUPROFEN},
        ]
    )
    assert summary.n_total == 2
    assert summary.n_valid == 1
    assert summary.n_invalid == 1
    assert records[1].lipinski_strict_pass is True


def test_properties_file_audit(tmp_path: Path):
    input_path = tmp_path / "input.csv"
    input_path.write_text(
        "compound_id,smiles,label\n"
        f"asp,{ASPIRIN},active\n"
        f"ste,{STEARIC},inactive\n",
        encoding="utf-8",
    )
    records, summary, manifest = run_properties_file(
        input_path, tmp_path / "out", PropertiesConfig()
    )
    assert (tmp_path / "out" / "properties.csv").is_file()
    audit = json.loads((tmp_path / "out" / "audit_log.json").read_text(encoding="utf-8"))
    assert audit["command"] == "properties"
    assert audit["hitpilot_version"] == "0.1.0"
    assert audit["rdkit_version"]
    assert audit["input_sha256"] == manifest.input_sha256
    assert audit["summary"]["n_total"] == 2
    assert summary.by_label["active"].lipinski_strict_pass_count == 1
    assert summary.by_label["inactive"].veber_pass_count == 0
    assert records[0].compound_id == "asp"
