"""Unit tests for Murcko scaffolds, Butina clustering, and sphere exclusion."""

from __future__ import annotations

import json
from pathlib import Path

from rdkit import Chem, DataStructs

from hitpilot.diversity import (
    murcko_scaffolds,
    run_diversity,
    run_diversity_file,
)
from hitpilot.fingerprints import fingerprint, morgan_generator
from hitpilot.models import DiversityConfig

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
IBUPROFEN = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
BENZENE = "c1ccccc1"
BENZENE_KEKULE = "C1=CC=CC=C1"
ETHANOL = "CCO"
HEXANE = "CCCCCC"
ERLOTINIB = "C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1"


def test_acyclic_molecules_have_empty_scaffolds():
    for smiles in (ETHANOL, HEXANE):
        murcko, generic, acyclic = murcko_scaffolds(Chem.MolFromSmiles(smiles))
        assert acyclic is True
        assert murcko == ""
        assert generic == ""


def test_murcko_and_generic_scaffolds_for_known_drugs():
    benzene_m, benzene_g, acyclic = murcko_scaffolds(Chem.MolFromSmiles(BENZENE))
    kekule_m, kekule_g, kekule_ac = murcko_scaffolds(Chem.MolFromSmiles(BENZENE_KEKULE))
    assert acyclic is False and kekule_ac is False
    assert benzene_m == kekule_m == "c1ccccc1"
    assert benzene_g == kekule_g == "C1CCCCC1"

    ibu_m, ibu_g, _ = murcko_scaffolds(Chem.MolFromSmiles(IBUPROFEN))
    asp_m, _, _ = murcko_scaffolds(Chem.MolFromSmiles(ASPIRIN))
    assert ibu_m == "c1ccccc1"
    assert asp_m == "c1ccccc1"
    assert ibu_g == "C1CCCCC1"

    erl_m, erl_g, erl_ac = murcko_scaffolds(Chem.MolFromSmiles(ERLOTINIB))
    assert erl_ac is False
    assert erl_m == "c1ccc(Nc2ncnc3ccccc23)cc1"
    assert erl_g == "C1CCC(CC2CCCC3CCCCC32)CC1"


def test_butina_clusters_and_scaffold_frequency():
    rows = [
        {"compound_id": "benzene", "smiles": BENZENE},
        {"compound_id": "kekule", "smiles": BENZENE_KEKULE},
        {"compound_id": "ibuprofen", "smiles": IBUPROFEN},
        {"compound_id": "ethanol", "smiles": ETHANOL},
    ]
    records, summary = run_diversity(rows, DiversityConfig(tanimoto_cutoff=0.60))
    recs = {rec.compound_id: rec for rec in records}
    # Identical benzene graphs share a cluster.
    assert recs["benzene"].cluster_id == recs["kekule"].cluster_id
    assert recs["benzene"].cluster_size == 2
    assert recs["benzene"].is_cluster_centroid or recs["kekule"].is_cluster_centroid
    assert recs["ibuprofen"].cluster_id != recs["benzene"].cluster_id
    assert recs["ethanol"].acyclic is True
    assert recs["benzene"].scaffold_frequency == recs["ibuprofen"].scaffold_frequency == 3
    assert summary.n_unique_murcko_scaffolds == 1  # benzene scaffold; ethanol is acyclic
    assert summary.n_acyclic == 1
    assert summary.n_clusters == 3


def test_sphere_exclusion_invariant_and_exclusion_links():
    rows = [
        {"compound_id": "benzene", "smiles": BENZENE, "label": "inactive"},
        {"compound_id": "kekule", "smiles": BENZENE_KEKULE, "label": "inactive"},
        {"compound_id": "ibuprofen", "smiles": IBUPROFEN, "label": "active"},
        {"compound_id": "aspirin", "smiles": ASPIRIN, "label": "active"},
        {"compound_id": "ethanol", "smiles": ETHANOL, "label": "inactive"},
    ]
    records, summary = run_diversity(
        rows, DiversityConfig(exclusion_similarity=0.65, tanimoto_cutoff=0.60)
    )
    selected = [rec for rec in records if rec.sphere_selected]
    assert summary.n_sphere_selected == len(selected)
    # Benzene and its kekule form cannot both be selected.
    ids = {rec.compound_id for rec in selected}
    assert "benzene" in ids  # first in input order is the leader
    assert "kekule" not in ids
    kekule = next(rec for rec in records if rec.compound_id == "kekule")
    assert kekule.excluded_by == "benzene"
    assert kekule.sphere_rank is None

    gen = morgan_generator()
    fps = {rec.compound_id: fingerprint(Chem.MolFromSmiles(rec.smiles), gen) for rec in selected}
    selected_ids = [rec.compound_id for rec in selected]
    for i, left in enumerate(selected_ids):
        for right in selected_ids[i + 1 :]:
            sim = DataStructs.TanimotoSimilarity(fps[left], fps[right])
            assert sim < 0.65, (left, right, sim)
    assert summary.n_actives_total == 2
    assert summary.n_actives_selected >= 1


def test_score_ranked_sphere_exclusion_keeps_best_leader():
    rows = [
        {"compound_id": "low", "smiles": BENZENE, "score": 0.2},
        {"compound_id": "high", "smiles": BENZENE_KEKULE, "score": 9.9},
        {"compound_id": "ibu", "smiles": IBUPROFEN, "score": 1.0},
    ]
    records, _ = run_diversity(
        rows,
        DiversityConfig(
            exclusion_similarity=0.65,
            rank_by_score=True,
            score_column="score",
            score_ascending=False,
        ),
    )
    recs = {rec.compound_id: rec for rec in records}
    assert recs["high"].sphere_selected is True
    assert recs["high"].sphere_rank == 1
    assert recs["low"].sphere_selected is False
    assert recs["low"].excluded_by == "high"


def test_diversity_file_outputs(tmp_path: Path):
    input_path = tmp_path / "in.csv"
    input_path.write_text(
        "compound_id,smiles,label\n"
        f"a,{BENZENE},inactive\n"
        f"b,{IBUPROFEN},active\n"
        f"c,{ETHANOL},inactive\n",
        encoding="utf-8",
    )
    records, summary, manifest = run_diversity_file(input_path, tmp_path / "out")
    assert (tmp_path / "out" / "diversity.csv").is_file()
    assert (tmp_path / "out" / "diverse_subset.csv").is_file()
    assert (tmp_path / "out" / "scaffold_summary.csv").is_file()
    audit = json.loads((tmp_path / "out" / "audit_log.json").read_text(encoding="utf-8"))
    assert audit["command"] == "diversity"
    assert audit["summary"]["n_total"] == 3
    assert summary.n_sphere_selected >= 1
    assert manifest.input_sha256
    assert any(rec.sphere_selected for rec in records)
