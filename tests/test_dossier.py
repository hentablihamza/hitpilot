"""Tests for the med-chem dossier, figures, and offline mini-benchmark."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from PIL import Image
from typer.testing import CliRunner

from hitpilot.benchmark import run_litpcba_benchmark
from hitpilot.cli import app
from hitpilot.dossier import (
    FIGURE_NAMES,
    DossierConfig,
    generate_medchem_dossier,
    generate_medchem_dossier_from_dirs,
)
from hitpilot.io_utils import sha256_file
from hitpilot.novelty import canonicalize_smiles

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
IBUPROFEN = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
BENZENE = "c1ccccc1"
ETHANOL = "CCO"
HEXANE = "CCCCCC"
CATECHOL = "Oc1ccccc1O"
NITRO = "O=[N+]([O-])c1ccccc1"
runner = CliRunner()


def _frame_shortlist() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "compound_id": ["ibu", "ben", "hex", "etoh"],
            "smiles": [IBUPROFEN, BENZENE, HEXANE, ETHANOL],
            "standardized_smiles": [IBUPROFEN, BENZENE, HEXANE, ETHANOL],
            "label": ["active", "inactive", "inactive", "inactive"],
            "score": [0.90, 0.40, 0.30, 0.10],
            "shortlist_rank": [1, 2, 3, 4],
            "mw": [206.285, 78.114, 86.178, 46.069],
            "clogp": [3.07, 1.69, 2.99, -0.16],
            "tpsa": [37.3, 0.0, 0.0, 20.23],
            "qed": [0.82, 0.45, 0.40, 0.39],
            "sa_score": [2.19, 1.00, 1.50, 1.20],
            "lipinski_pass": ["True", "True", "True", "True"],
            "veber_pass": ["True", "True", "True", "True"],
            "max_tanimoto": [1.0, 0.12, 0.08, 0.05],
            "nearest_ref_id": ["CHEMBL521", "CHEMBL1", "CHEMBL2", "CHEMBL3"],
            "nearest_ref_pchembl": [7.2, 6.1, 6.0, 6.0],
            "novelty_tier": [
                "known_or_close_match",
                "novel_chemotype",
                "novel_chemotype",
                "novel_chemotype",
            ],
            "murcko_scaffold_smiles": ["c1ccccc1", "c1ccccc1", "", ""],
            "cluster_id": ["1", "1", "2", "3"],
            "cluster_size": [2, 2, 1, 1],
        }
    )


def _frame_all() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "compound_id": ["ibu", "ben", "hex", "etoh", "asp", "cat", "nitro"],
            "smiles": [IBUPROFEN, BENZENE, HEXANE, ETHANOL, ASPIRIN, CATECHOL, NITRO],
            "standardized_smiles": [IBUPROFEN, BENZENE, HEXANE, ETHANOL, ASPIRIN, CATECHOL, NITRO],
            "inchikey": ["IBU", "BEN", "HEX", "ETO", "ASP", "CAT", "NIT"],
            "passed": ["True", "True", "True", "True", "False", "False", "False"],
            "label": [
                "active",
                "inactive",
                "inactive",
                "inactive",
                "inactive",
                "inactive",
                "active",
            ],
            "alerts": [
                "",
                "",
                "",
                "",
                "BRENK:phenol_ester",
                "PAINS_A:catechol_A(92)",
                "BRENK:nitro_group|NIH:nitro_group",
            ],
            "drop_reasons": ["", "", "", "", "brenk_alert", "pains_alert", "brenk_alert|nih_alert"],
        }
    )


def _frame_phase2(short: pd.DataFrame) -> pd.DataFrame:
    extra = pd.DataFrame(
        {
            "compound_id": ["extra"],
            "smiles": ["CCC"],
            "standardized_smiles": ["CCC"],
            "label": ["inactive"],
            "mw": [44.1],
            "clogp": [1.4],
            "tpsa": [0.0],
            "qed": [0.41],
            "sa_score": [1.1],
            "max_tanimoto": [0.04],
            "novelty_tier": ["novel_chemotype"],
            "murcko_scaffold_smiles": [""],
            "cluster_id": ["4"],
            "cluster_size": [1],
        }
    )
    return pd.concat([short, extra], ignore_index=True)


def _hygiene_audit() -> dict:
    return {
        "summary": {
            "n_total": 7,
            "n_passed": 4,
            "n_dropped": 3,
            "pass_rate": 4 / 7,
            "counts_by_drop_reason": {"brenk_alert": 2, "pains_alert": 1, "nih_alert": 1},
            "counts_by_catalog": {"BRENK": 2, "PAINS_A": 1, "NIH": 1},
            "by_label": {
                "active": {"n_total": 2, "n_passed": 1, "n_dropped": 1, "pass_rate": 0.5},
                "inactive": {"n_total": 5, "n_passed": 3, "n_dropped": 2, "pass_rate": 0.6},
            },
        }
    }


def _scoring_audit() -> dict:
    return {
        "n_matched": 4,
        "n_shortlist": 4,
        "n_shortlist_scaffolds": 2,
        "diverse_hits_top_n": 1,
        "higher_is_better": True,
        "enrichment": {
            "n_compounds": 4,
            "n_actives": 1,
            "n_inactives": 3,
            "roc_auc": 1.0,
            "pr_auc": 1.0,
            "ef_1pct": 4.0,
            "ef_5pct": 4.0,
            "ef_10pct": 4.0,
            "bedroc_20": 1.0,
            "bedroc_80_5": 1.0,
            "hits_top_25": 1,
            "hits_top_50": 1,
            "hits_top_100": 1,
            "k_1pct": 1,
            "k_5pct": 1,
            "k_10pct": 1,
            "hits_at_1pct": 1,
            "hits_at_5pct": 1,
            "hits_at_10pct": 1,
        },
    }


def test_dossier_writes_pngs_markdown_html_and_manifest(tmp_path: Path):
    short = _frame_shortlist()
    all_compounds = _frame_all()
    phase2 = _frame_phase2(short)
    scaf = pd.DataFrame(
        {
            "murcko_scaffold_smiles": ["c1ccccc1", ""],
            "generic_scaffold_smiles": ["C1CCCCC1", ""],
            "acyclic": [False, True],
            "count": [2, 2],
            "n_active": [1, 0],
            "n_inactive": [1, 2],
            "example_compound_ids": ["ibu|ben", "hex|etoh"],
        }
    )
    # Persist sources so the manifest can hash them.
    all_path = tmp_path / "all_compounds.csv"
    p2_path = tmp_path / "phase2_enriched.csv"
    ranked_path = tmp_path / "ranked.csv"
    short_path = tmp_path / "shortlist.csv"
    scaf_path = tmp_path / "scaffolds.csv"
    all_compounds.to_csv(all_path, index=False)
    phase2.to_csv(p2_path, index=False)
    short.to_csv(ranked_path, index=False)
    short.to_csv(short_path, index=False)
    scaf.to_csv(scaf_path, index=False)
    out = tmp_path / "dossier"
    result = generate_medchem_dossier(
        out,
        all_compounds=all_path,
        phase2_enriched=p2_path,
        ranked=ranked_path,
        shortlist=short_path,
        hygiene_audit=_hygiene_audit(),
        phase2_audit={
            "summary": {
                "n_input": 4,
                "n_reference_actives": 1,
                "properties": {},
                "novelty": {},
                "diversity": {},
            }
        },
        scoring_audit=_scoring_audit(),
        baseline_stats={
            "n_overlap_removed": 1,
            "train_standardized": {"n": 8, "n_active": 2, "n_inactive": 6},
            "valid_standardized": {"n": 5, "n_active": 2, "n_inactive": 3},
            "valid_leakfree": {"n": 4, "n_active": 1, "n_inactive": 3},
        },
        scaffold_summary=scaf_path,
        config=DossierConfig(target="ESR1_ant", target_chembl_id="CHEMBL206", top_n=4, n_grid=4),
    )
    assert Path(result.markdown_path).is_file()
    assert Path(result.html_path).is_file()
    assert Path(result.manifest_path).is_file()
    markdown = Path(result.markdown_path).read_text(encoding="utf-8")
    html = Path(result.html_path).read_text(encoding="utf-8")
    assert "Leak-free" in markdown
    assert "**Standardized SMILES overlap removed from `_V`:** 1" in markdown
    assert "ROC-AUC" in markdown
    assert "ibu" in markdown
    assert "figures/top_shortlist_grid.png" in markdown
    assert "data:image/png;base64," in html
    assert html.count("data:image/png;base64,") == 5
    assert "src='http" not in html and 'src="http' not in html
    assert "badge-active" in html
    for name in FIGURE_NAMES:
        path = Path(result.figures[name])
        assert path.is_file(), name
        data = path.read_bytes()
        assert data.startswith(b"\x89PNG\r\n\x1a\n"), name
        with Image.open(path) as image:
            assert image.format == "PNG"
            width, height = image.size
            assert width >= 400, (name, image.size)
            assert height >= 200, (name, image.size)
        # Grid of 4 columns should be at least close to 4 * 380 for shortlist.
        if name == "top_shortlist_grid.png":
            with Image.open(path) as image:
                assert image.size[0] >= 1400

    manifest = json.loads(Path(result.manifest_path).read_text(encoding="utf-8"))
    assert manifest["command"] == "dossier"
    assert manifest["summary"]["baseline"]["n_overlap_removed"] == 1
    assert manifest["summary"]["enrichment"]["roc_auc"] == 1.0
    for name, meta in manifest["artifacts"].items():
        artifact = Path(meta["path"])
        assert artifact.is_file()
        if name == "dossier_manifest.json":
            continue
        assert meta["sha256"] == sha256_file(artifact)
        assert meta["bytes"] == artifact.stat().st_size
    for name, meta in manifest["inputs"].items():
        assert meta["sha256"] == sha256_file(meta["path"])


def _write_mini_litpcba(root: Path) -> None:
    target = root / "ESR1_ant"
    raw = target / "raw"
    raw.mkdir(parents=True)
    (raw / "AID_active_T.smi").write_text(f"{IBUPROFEN} T-IBU\nCC(=O)O T-HOAC\n", encoding="utf-8")
    (raw / "AID_inactive_T.smi").write_text(
        f"{ASPIRIN} T-ASP\n{HEXANE} T-HEX\n{BENZENE} T-BENZ\n{ETHANOL} T-ETOH\n"
        f"{CATECHOL} T-CAT\nCCC T-PROP\n",
        encoding="utf-8",
    )
    (raw / "AID_active_V.smi").write_text(
        f"c1ccc(CC(=O)O)cc1 V-PAA\n{IBUPROFEN} V-IBU-LEAK\n",
        encoding="utf-8",
    )
    (raw / "AID_inactive_V.smi").write_text(
        "CCCC V-BUT\nCCOCC V-ETHER\n",
        encoding="utf-8",
    )
    rows = [
        ("T-IBU", IBUPROFEN, "active"),
        ("T-HOAC", "CC(=O)O", "active"),
        ("T-ASP", ASPIRIN, "inactive"),
        ("T-HEX", HEXANE, "inactive"),
        ("T-BENZ", BENZENE, "inactive"),
        ("T-ETOH", ETHANOL, "inactive"),
        ("T-CAT", CATECHOL, "inactive"),
        ("T-PROP", "CCC", "inactive"),
        ("V-PAA", "c1ccc(CC(=O)O)cc1", "active"),
        ("V-IBU-LEAK", IBUPROFEN, "active"),
        ("V-BUT", "CCCC", "inactive"),
        ("V-ETHER", "CCOCC", "inactive"),
    ]
    csv = target / "esr1_ant_litpcba.csv"
    csv.write_text(
        "compound_id,smiles,label\n" + "".join(f"{cid},{smi},{lab}\n" for cid, smi, lab in rows),
        encoding="utf-8",
    )


def _write_chembl_cache(cache_dir: Path, target_chembl_id: str = "CHEMBL206") -> None:
    cache_dir.mkdir(parents=True, exist_ok=True)
    canon = canonicalize_smiles(IBUPROFEN)
    pd.DataFrame(
        [
            {
                "molecule_chembl_id": "CHEMBL521",
                "smiles": IBUPROFEN,
                "canonical_smiles": canon,
                "pchembl_value": 7.2,
                "target_chembl_id": target_chembl_id,
            }
        ]
    ).to_csv(cache_dir / f"{target_chembl_id}_actives.csv", index=False)
    (cache_dir / f"{target_chembl_id}_actives.meta.json").write_text(
        json.dumps(
            {
                "target_chembl_id": target_chembl_id,
                "pchembl_min": 6.0,
                "max_records": 2000,
                "n_activity_records": 1,
                "n_unique_actives": 1,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


def test_offline_mini_benchmark_and_cli(tmp_path: Path):
    data_dir = tmp_path / "litpcba"
    cache_dir = tmp_path / "chembl"
    _write_mini_litpcba(data_dir)
    _write_chembl_cache(cache_dir)
    out = tmp_path / "benchmark"
    result = run_litpcba_benchmark(
        target="ESR1_ant",
        output_dir=out,
        target_chembl_id="CHEMBL206",
        data_dir=data_dir,
        chembl_cache_dir=cache_dir,
        top_n=3,
        n_estimators=20,
        random_state=42,
    )
    assert result["cached_litpcba"] is True
    assert (out / "phase1" / "passed_compounds.csv").is_file()
    assert (out / "phase2" / "phase2_enriched.csv").is_file()
    assert (out / "phase3" / "baseline_ecfp4_rf_scores.csv").is_file()
    assert (out / "phase3" / "baseline_stats.json").is_file()
    assert (out / "dossier" / "medchem_dossier.md").is_file()
    assert (out / "dossier" / "medchem_dossier.html").is_file()
    assert (out / "dossier" / "dossier_manifest.json").is_file()
    assert result["baseline"]["n_overlap_removed"] >= 1
    html = (out / "dossier" / "medchem_dossier.html").read_text(encoding="utf-8")
    assert html.count("data:image/png;base64,") == 5
    for name in FIGURE_NAMES:
        path = out / "dossier" / "figures" / name
        assert path.is_file()
        assert path.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")

    # CLI dossier command on the same artifacts.
    dossier_out = tmp_path / "dossier_cli"
    cli = runner.invoke(
        app,
        [
            "dossier",
            "--phase1-dir",
            str(out / "phase1"),
            "--phase2-dir",
            str(out / "phase2"),
            "--phase3-dir",
            str(out / "phase3"),
            "--output-dir",
            str(dossier_out),
            "--target",
            "ESR1_ant",
            "--top-n",
            "3",
        ],
        catch_exceptions=False,
    )
    assert cli.exit_code == 0, cli.output
    assert (dossier_out / "medchem_dossier.md").is_file()

    # CLI benchmark on a fresh output dir (still cached / offline).
    bench_out = tmp_path / "benchmark_cli"
    cli = runner.invoke(
        app,
        [
            "benchmark",
            "ESR1_ant",
            "--output-dir",
            str(bench_out),
            "--data-dir",
            str(data_dir),
            "--chembl-cache-dir",
            str(cache_dir),
            "--top-n",
            "3",
        ],
        catch_exceptions=False,
    )
    assert cli.exit_code == 0, cli.output
    assert (bench_out / "dossier" / "medchem_dossier.html").is_file()


def test_generate_from_dirs_roundtrip(tmp_path: Path):
    # Minimal dirs so from_dirs can be unit-tested without a full benchmark.
    p1 = tmp_path / "phase1"
    p2 = tmp_path / "phase2"
    p3 = tmp_path / "phase3"
    for folder in (p1, p2, p3):
        folder.mkdir()
    _frame_all().to_csv(p1 / "all_compounds.csv", index=False)
    (p1 / "audit_log.json").write_text(json.dumps(_hygiene_audit()), encoding="utf-8")
    short = _frame_shortlist()
    _frame_phase2(short).to_csv(p2 / "phase2_enriched.csv", index=False)
    (p2 / "audit_log.json").write_text(json.dumps({"summary": {"n_input": 4}}), encoding="utf-8")
    pd.DataFrame(
        {
            "murcko_scaffold_smiles": ["c1ccccc1"],
            "generic_scaffold_smiles": ["C1CCCCC1"],
            "acyclic": [False],
            "count": [2],
            "n_active": [1],
            "n_inactive": [1],
            "example_compound_ids": ["ibu"],
        }
    ).to_csv(p2 / "scaffold_summary.csv", index=False)
    short.to_csv(p3 / "ranked_compounds.csv", index=False)
    short.to_csv(p3 / "shortlist_topN_diverse.csv", index=False)
    (p3 / "audit_log.json").write_text(json.dumps(_scoring_audit()), encoding="utf-8")
    (p3 / "baseline_stats.json").write_text(
        json.dumps(
            {"n_overlap_removed": 2, "valid_leakfree": {"n": 4, "n_active": 1, "n_inactive": 3}}
        ),
        encoding="utf-8",
    )
    result = generate_medchem_dossier_from_dirs(
        tmp_path / "out",
        phase1_dir=p1,
        phase2_dir=p2,
        phase3_dir=p3,
        config=DossierConfig(top_n=4, n_grid=4),
    )
    assert Path(result.html_path).is_file()
    assert result.manifest.summary["baseline"]["n_overlap_removed"] == 2
