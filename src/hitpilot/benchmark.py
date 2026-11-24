"""End-to-end LIT-PCBA public benchmark orchestrator (Phases 1-4)."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from hitpilot.data.litpcba import (
    KNOWN_TARGETS,
    fetch_litpcba,
    load_litpcba,
    normalize_target,
    unified_csv_path,
    write_unified_csv,
)
from hitpilot.dossier import DossierConfig, generate_medchem_dossier_from_dirs
from hitpilot.hygiene import run_hygiene_file
from hitpilot.io_utils import write_audit
from hitpilot.models import DiversityConfig, HygieneConfig, NoveltyConfig, PropertiesConfig
from hitpilot.phase2 import run_phase2_file
from hitpilot.scoring import (
    ScorerHookConfig,
    generate_litpcba_baseline_scores,
    run_scorer_hook_file,
)


def _jsonable(value: Any) -> Any:
    if hasattr(value, "model_dump"):
        return value.model_dump(mode="json")
    if isinstance(value, dict):
        return {key: _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, pd.DataFrame):
        return {"n_rows": int(len(value)), "columns": list(value.columns)}
    return value


def load_or_fetch_litpcba(
    target: str,
    data_dir: str | Path = "data/litpcba",
) -> dict[str, Any]:
    """Reuse a cached LIT-PCBA target directory, fetching only if necessary."""
    target_name = normalize_target(target)
    data_dir = Path(data_dir)
    target_dir = data_dir / target_name
    csv_path = unified_csv_path(target_dir, target_name)
    raw_dir = target_dir / "raw"
    has_raw = raw_dir.is_dir() and any(raw_dir.rglob("*.smi"))
    if csv_path.is_file() or has_raw:
        frame = load_litpcba(target_dir, target=target_name)
        if not csv_path.is_file():
            csv_path = write_unified_csv(frame, csv_path)
        n_active = int((frame["label"] == "active").sum()) if "label" in frame.columns else 0
        n_inactive = int((frame["label"] == "inactive").sum()) if "label" in frame.columns else 0
        return {
            "target": target_name,
            "cached": True,
            "url": None,
            "target_dir": str(target_dir.resolve()),
            "csv_path": str(csv_path.resolve()),
            "n_compounds": int(len(frame)),
            "n_active": n_active,
            "n_inactive": n_inactive,
            "protein": KNOWN_TARGETS.get(target_name, {}).get("protein"),
            "description": KNOWN_TARGETS.get(target_name, {}).get("description"),
        }
    meta = fetch_litpcba(target_name, dest=data_dir)
    meta["cached"] = False
    return meta


def run_litpcba_benchmark(
    target: str = "ESR1_ant",
    output_dir: str | Path = "outputs/benchmark_ESR1_ant",
    *,
    target_chembl_id: str = "CHEMBL206",
    data_dir: str | Path = "data/litpcba",
    chembl_cache_dir: str | Path = "data/chembl",
    scores_csv: str | Path | None = None,
    top_n: int = 25,
    exclusion_similarity: float = 0.65,
    pchembl_min: float = 6.0,
    max_records: int = 2000,
    hygiene_config: HygieneConfig | None = None,
    properties_config: PropertiesConfig | None = None,
    novelty_config: NoveltyConfig | None = None,
    diversity_config: DiversityConfig | None = None,
    scorer_config: ScorerHookConfig | None = None,
    random_state: int = 42,
    n_estimators: int = 200,
    refresh_chembl_cache: bool = False,
) -> dict[str, Any]:
    """Run Phase 1 hygiene through Phase 4 dossier for one public LIT-PCBA target."""
    target_name = normalize_target(target)
    output_dir = Path(output_dir)
    phase1_dir = output_dir / "phase1"
    phase2_dir = output_dir / "phase2"
    phase3_dir = output_dir / "phase3"
    dossier_dir = output_dir / "dossier"
    for folder in (phase1_dir, phase2_dir, phase3_dir, dossier_dir):
        folder.mkdir(parents=True, exist_ok=True)

    litpcba = load_or_fetch_litpcba(target_name, data_dir=data_dir)
    input_csv = Path(litpcba["csv_path"])
    hygiene_config = hygiene_config or HygieneConfig()
    hygiene_manifest = run_hygiene_file(input_csv, phase1_dir, hygiene_config)

    novelty_config = novelty_config or NoveltyConfig(
        target_chembl_id=target_chembl_id,
        pchembl_min=pchembl_min,
        max_records=max_records,
        cache_dir=str(chembl_cache_dir),
        refresh_cache=refresh_chembl_cache,
    )
    properties_config = properties_config or PropertiesConfig()
    diversity_config = diversity_config or DiversityConfig(
        exclusion_similarity=exclusion_similarity,
    )
    phase2_result = run_phase2_file(
        phase1_dir / "passed_compounds.csv",
        phase2_dir,
        properties_config=properties_config,
        novelty_config=novelty_config,
        diversity_config=diversity_config,
    )

    scorer_config = scorer_config or ScorerHookConfig(
        top_n=top_n,
        diversity_cutoff=exclusion_similarity,
    )
    baseline_stats: dict[str, Any] | None = None
    if scores_csv is None:
        scores_path = phase3_dir / "baseline_ecfp4_rf_scores.csv"
        baseline_stats = generate_litpcba_baseline_scores(
            target=target_name,
            data_dir=data_dir,
            output_path=scores_path,
            random_state=random_state,
            n_estimators=n_estimators,
        )
        # Drop the large scores DataFrame before writing JSON.
        baseline_for_disk = {
            key: _jsonable(value) for key, value in baseline_stats.items() if key != "scores"
        }
        write_audit(phase3_dir / "baseline_stats.json", baseline_for_disk)
        baseline_stats = baseline_for_disk
    else:
        scores_path = Path(scores_csv)

    ranked, shortlist, scoring_summary, scoring_manifest = run_scorer_hook_file(
        phase2_dir / "phase2_enriched.csv",
        scores_path,
        phase3_dir,
        scorer_config,
    )

    description = (
        litpcba.get("description") or KNOWN_TARGETS.get(target_name, {}).get("description") or ""
    )
    dossier = generate_medchem_dossier_from_dirs(
        dossier_dir,
        phase1_dir=phase1_dir,
        phase2_dir=phase2_dir,
        phase3_dir=phase3_dir,
        config=DossierConfig(
            target=target_name,
            target_chembl_id=target_chembl_id,
            target_description=str(description),
            top_n=top_n,
        ),
    )

    result = {
        "target": target_name,
        "target_chembl_id": target_chembl_id,
        "output_dir": str(output_dir.resolve()),
        "litpcba": {key: _jsonable(value) for key, value in litpcba.items()},
        "phase1_dir": str(phase1_dir.resolve()),
        "phase2_dir": str(phase2_dir.resolve()),
        "phase3_dir": str(phase3_dir.resolve()),
        "dossier_dir": str(dossier_dir.resolve()),
        "hygiene": hygiene_manifest.model_dump(mode="json")
        if hasattr(hygiene_manifest, "model_dump")
        else hygiene_manifest,
        "phase2": _jsonable(phase2_result["manifest"].summary)
        if "manifest" in phase2_result
        else {},
        "baseline": baseline_stats,
        "scoring": scoring_summary.model_dump(mode="json"),
        "dossier": dossier.model_dump(mode="json"),
        "n_ranked": int(len(ranked)),
        "n_shortlist": int(len(shortlist)),
        "cached_litpcba": bool(litpcba.get("cached")),
    }
    write_audit(output_dir / "benchmark_manifest.json", result)
    return result
