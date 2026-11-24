"""Unit tests for the pluggable scorer hook, enrichment metrics, and baseline."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest
from rdkit import Chem, DataStructs
from rdkit.ML.Scoring import Scoring

from hitpilot.fingerprints import fingerprint, morgan_generator
from hitpilot.scoring import (
    ScorerHookConfig,
    compute_enrichment_metrics,
    generate_litpcba_baseline_scores,
    merge_external_scores,
    run_scorer_hook,
    run_scorer_hook_file,
)

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
IBUPROFEN = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
BENZENE = "c1ccccc1"
ETHANOL = "CCO"
HEXANE = "CCCCCC"


def _perfect_df(n: int = 100, n_actives: int = 5, higher_is_better: bool = True) -> pd.DataFrame:
    rows = []
    for i in range(n):
        is_active = i < n_actives
        if higher_is_better:
            score = 1.0 - i / n
        else:
            score = -10.0 + i / n  # more negative is better
        rows.append(
            {
                "compound_id": f"c{i}",
                "score": score,
                "label": "active" if is_active else "inactive",
            }
        )
    return pd.DataFrame(rows)


def test_perfect_ranking_metrics_higher_is_better():
    metrics = compute_enrichment_metrics(
        _perfect_df(100, 5, True), "score", higher_is_better=True
    )
    assert metrics is not None
    assert metrics.n_compounds == 100
    assert metrics.n_actives == 5
    assert metrics.n_inactives == 95
    assert metrics.roc_auc == pytest.approx(1.0)
    assert metrics.pr_auc == pytest.approx(1.0)
    assert metrics.k_1pct == 1
    assert metrics.hits_at_1pct == 1
    assert metrics.ef_1pct == pytest.approx(20.0)
    assert metrics.k_5pct == 5
    assert metrics.hits_at_5pct == 5
    assert metrics.ef_5pct == pytest.approx(20.0)
    assert metrics.k_10pct == 10
    assert metrics.hits_at_10pct == 5
    assert metrics.ef_10pct == pytest.approx(10.0)
    assert metrics.hits_top_25 == 5
    assert metrics.hits_top_50 == 5
    assert metrics.hits_top_100 == 5
    assert metrics.bedroc_20 == pytest.approx(1.0)
    assert metrics.bedroc_80_5 == pytest.approx(1.0)


def test_perfect_ranking_metrics_docking_lower_is_better():
    metrics = compute_enrichment_metrics(
        _perfect_df(100, 5, False), "score", higher_is_better=False
    )
    assert metrics is not None
    assert metrics.roc_auc == pytest.approx(1.0)
    assert metrics.pr_auc == pytest.approx(1.0)
    assert metrics.bedroc_20 == pytest.approx(1.0)
    assert metrics.bedroc_80_5 == pytest.approx(1.0)
    assert metrics.ef_1pct == pytest.approx(20.0)
    # Inverting the direction must collapse early recognition.
    inverted = compute_enrichment_metrics(
        _perfect_df(100, 5, False), "score", higher_is_better=True
    )
    assert inverted is not None
    assert inverted.roc_auc == pytest.approx(0.0)
    assert inverted.hits_top_25 == 0
    assert inverted.bedroc_20 < 0.1


def test_bedroc_matches_rdkit_official_implementation():
    labels = [[1], [1], [0], [0], [0], [0], [0], [0], [0], [0]]
    assert Scoring.CalcBEDROC(labels, 0, 20.0) == pytest.approx(1.0)
    df = pd.DataFrame(
        {
            "score": [10, 9, 1, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2],
            "label": ["active", "active"] + ["inactive"] * 8,
        }
    )
    metrics = compute_enrichment_metrics(df, "score", higher_is_better=True)
    assert metrics is not None
    assert metrics.bedroc_20 == pytest.approx(Scoring.CalcBEDROC(labels, 0, 20.0))
    assert metrics.bedroc_80_5 == pytest.approx(Scoring.CalcBEDROC(labels, 0, 80.5))


def test_enrichment_requires_both_classes():
    df = pd.DataFrame({"score": [1.0, 0.5], "label": ["active", "active"]})
    assert compute_enrichment_metrics(df, "score") is None


def test_merge_by_id_and_canonical_smiles_fallback():
    compounds = pd.DataFrame(
        {
            "compound_id": ["A", "B", "C"],
            "smiles": [IBUPROFEN, ASPIRIN, ETHANOL],
            "label": ["active", "inactive", "inactive"],
        }
    )
    scores = pd.DataFrame(
        {
            "compound_id": ["A", "SMILES-ONLY", "UNMATCHED-SCORE"],
            "smiles": [IBUPROFEN, "CC(=O)Oc1ccccc1C(=O)O", HEXANE],
            "score": [0.9, 0.4, 0.1],
        }
    )
    ranked, stats = merge_external_scores(compounds, scores, ScorerHookConfig())
    assert stats["n_matched_by_id"] == 1
    assert stats["n_matched_by_smiles"] == 1
    assert stats["n_matched"] == 2
    assert stats["n_unmatched_compounds"] == 1  # ethanol
    assert stats["n_unmatched_scores"] == 1  # hexane / UNMATCHED-SCORE
    recs = {row.compound_id: row for row in ranked.itertuples()}
    assert recs["A"].score_rank == 1
    assert recs["B"].score_rank == 2
    assert "C" not in recs
    # higher score first
    assert list(ranked["compound_id"]) == ["A", "B"]


def test_merge_lower_is_better_ranks_most_negative_first():
    compounds = pd.DataFrame(
        {
            "compound_id": ["good", "bad"],
            "smiles": [IBUPROFEN, ETHANOL],
            "label": ["active", "inactive"],
        }
    )
    scores = pd.DataFrame({"compound_id": ["good", "bad"], "score": [-12.4, -3.1]})
    ranked, _stats = merge_external_scores(
        compounds, scores, ScorerHookConfig(higher_is_better=False)
    )
    assert list(ranked["compound_id"]) == ["good", "bad"]
    assert list(ranked["score_rank"]) == [1, 2]


def test_score_ordered_sphere_exclusion_invariant():
    compounds = pd.DataFrame(
        {
            "compound_id": ["benzene", "kekule", "ibuprofen", "ethanol", "hexane"],
            "smiles": [BENZENE, "C1=CC=CC=C1", IBUPROFEN, ETHANOL, HEXANE],
            "score": [0.99, 0.98, 0.80, 0.50, 0.40],
            "label": ["inactive", "inactive", "active", "inactive", "inactive"],
        }
    )
    scores = compounds[["compound_id", "score"]]
    ranked, shortlist, summary, _records = run_scorer_hook(
        compounds,
        scores,
        ScorerHookConfig(top_n=3, diversity_cutoff=0.65, higher_is_better=True),
    )
    assert list(ranked["compound_id"]) == [
        "benzene",
        "kekule",
        "ibuprofen",
        "ethanol",
        "hexane",
    ]
    ids = list(shortlist["compound_id"])
    assert "benzene" in ids
    assert "kekule" not in ids  # same structure, lower score, excluded
    assert "ibuprofen" in ids
    assert len(shortlist) == 3
    assert summary.diverse_hits_top_n == 1
    gen = morgan_generator()
    fps = {
        rec.compound_id: fingerprint(Chem.MolFromSmiles(rec.smiles), gen)
        for rec in shortlist.itertuples()
    }
    chosen = list(shortlist["compound_id"])
    for i, left in enumerate(chosen):
        for right in chosen[i + 1 :]:
            sim = DataStructs.TanimotoSimilarity(fps[left], fps[right])
            assert sim < 0.65, (left, right, sim)


def test_scorer_hook_file_audit_checksums_both_inputs(tmp_path: Path):
    compounds = tmp_path / "compounds.csv"
    scores = tmp_path / "scores.csv"
    compounds.write_text(
        "compound_id,smiles,label\n"
        f"ibu,{IBUPROFEN},active\n"
        f"asp,{ASPIRIN},inactive\n",
        encoding="utf-8",
    )
    scores.write_text("compound_id,score,model_name,split\nibu,0.8,ecfp4_rf,V\nasp,0.2,ecfp4_rf,V\n")
    ranked, shortlist, summary, manifest = run_scorer_hook_file(
        compounds, scores, tmp_path / "out", ScorerHookConfig(top_n=2)
    )
    assert (tmp_path / "out" / "ranked_compounds.csv").is_file()
    assert (tmp_path / "out" / "shortlist_topN_diverse.csv").is_file()
    audit = json.loads((tmp_path / "out" / "audit_log.json").read_text(encoding="utf-8"))
    assert audit["command"] == "score"
    assert audit["config"]["compounds_sha256"] == manifest.input_sha256
    assert audit["config"]["scores_sha256"]
    assert audit["config"]["scores_sha256"] != audit["config"]["compounds_sha256"]
    assert audit["summary"]["enrichment"]["roc_auc"] == pytest.approx(1.0)
    assert summary.n_matched == 2
    assert len(shortlist) == 2
    assert list(ranked["compound_id"]) == ["ibu", "asp"]


def _write_smi(path: Path, rows: list[tuple[str, str]]) -> None:
    path.write_text("".join(f"{smi} {cid}\n" for smi, cid in rows), encoding="utf-8")


def test_leak_free_baseline_removes_tv_smiles_overlap(tmp_path: Path):
    # Distinct actives in T vs V, plus a leaked ethanol that appears in both.
    root = tmp_path / "ESR1_ant"
    root.mkdir()
    _write_smi(root / "AID_active_T.smi", [(IBUPROFEN, "T-IBU"), (ETHANOL, "T-ETOH")])
    _write_smi(root / "AID_inactive_T.smi", [(ASPIRIN, "T-ASP"), (HEXANE, "T-HEX")])
    _write_smi(
        root / "AID_active_V.smi",
        [(BENZENE, "V-BENZ"), ("CCO", "V-ETOH-LEAK")],  # CCO == ethanol
    )
    _write_smi(root / "AID_inactive_V.smi", [("CCC", "V-PROP")])
    output = tmp_path / "scores.csv"
    result = generate_litpcba_baseline_scores(
        target="ESR1_ant",
        data_dir=tmp_path,
        output_path=output,
        random_state=42,
        n_estimators=20,
    )
    assert result["n_overlap_removed"] == 1
    assert "V-ETOH-LEAK" in result["overlap_ids"]
    scores = pd.read_csv(output, dtype={"compound_id": str})
    assert set(scores.columns) >= {"compound_id", "score", "model_name", "split"}
    assert set(scores["compound_id"]) == {"V-BENZ", "V-PROP"}
    assert "V-ETOH-LEAK" not in set(scores["compound_id"])
    assert set(scores["split"]) == {"V"}
    assert set(scores["model_name"]) == {"ecfp4_rf"}
    assert result["valid_leakfree"]["n"] == 2
    assert result["train_standardized"]["n"] == 4
    assert scores["score"].between(0.0, 1.0).all()



def test_merge_overrides_existing_empty_score_column():
    compounds = pd.DataFrame(
        {
            "compound_id": ["A", "B"],
            "smiles": [IBUPROFEN, ETHANOL],
            "label": ["active", "inactive"],
            "score": ["", ""],
        }
    )
    scores = pd.DataFrame({"compound_id": ["A", "B"], "score": [0.9, 0.1]})
    ranked, stats = merge_external_scores(compounds, scores, ScorerHookConfig())
    assert stats["n_matched"] == 2
    assert list(ranked["compound_id"]) == ["A", "B"]
    assert list(ranked["score"]) == [0.9, 0.1]
    metrics = compute_enrichment_metrics(ranked, "score", higher_is_better=True)
    assert metrics is not None
    assert metrics.roc_auc == pytest.approx(1.0)
