"""Orchestrate Phase 2 properties, novelty, and diversity into one audit trail."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from hitpilot.diversity import run_diversity, scaffold_summary_rows
from hitpilot.io_utils import build_audit, load_compound_frame, write_audit, write_csv
from hitpilot.models import DiversityConfig, NoveltyConfig, PropertiesConfig
from hitpilot.novelty import fetch_chembl_target_actives, run_novelty
from hitpilot.properties import run_properties


def _drop_overlap(
    frame: pd.DataFrame,
    new_cols: list[str],
    keep: tuple[str, ...] = ("compound_id",),
) -> pd.DataFrame:
    drop = [col for col in new_cols if col in frame.columns and col not in keep]
    return frame.drop(columns=drop, errors="ignore")


def run_phase2_file(
    input_path: str | Path,
    output_dir: str | Path,
    *,
    properties_config: PropertiesConfig | None = None,
    novelty_config: NoveltyConfig | None = None,
    diversity_config: DiversityConfig | None = None,
    smiles_column: str | None = None,
    id_column: str | None = None,
    label_column: str | None = None,
    score_column: str | None = None,
):
    """Run properties + novelty + diversity and write unified Phase 2 outputs."""
    properties_config = properties_config or PropertiesConfig()
    novelty_config = novelty_config or NoveltyConfig()
    diversity_config = diversity_config or DiversityConfig()
    if score_column and not diversity_config.score_column:
        diversity_config = diversity_config.model_copy(
            update={"score_column": score_column, "rank_by_score": True}
        )

    input_path = Path(input_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    frame = load_compound_frame(
        input_path,
        smiles_column=smiles_column,
        id_column=id_column,
        label_column=label_column,
        score_column=score_column or diversity_config.score_column,
    )
    rows = frame.to_dict(orient="records")

    references = fetch_chembl_target_actives(
        target_chembl_id=novelty_config.target_chembl_id,
        pchembl_min=novelty_config.pchembl_min,
        max_records=novelty_config.max_records,
        cache_dir=novelty_config.cache_dir,
        refresh_cache=novelty_config.refresh_cache,
        timeout=novelty_config.timeout,
    )
    property_records, property_summary = run_properties(rows, properties_config)
    novelty_records, novelty_summary = run_novelty(rows, references, novelty_config)
    diversity_records, diversity_summary = run_diversity(rows, diversity_config)

    prop_df = pd.DataFrame([rec.to_csv_row() for rec in property_records])
    nov_df = pd.DataFrame([rec.to_csv_row() for rec in novelty_records])
    div_df = pd.DataFrame([rec.to_csv_row() for rec in diversity_records])

    enriched = frame.copy()
    enriched = _drop_overlap(enriched, [c for c in prop_df.columns if c != "compound_id"])
    enriched = enriched.merge(prop_df, on="compound_id", how="left", suffixes=("", "_prop"))
    nov_keep = {"compound_id", "smiles", "label"}
    enriched = _drop_overlap(enriched, [c for c in nov_df.columns if c not in nov_keep])
    enriched = enriched.merge(
        nov_df.drop(columns=[c for c in ("smiles", "label") if c in nov_df.columns]),
        on="compound_id",
        how="left",
    )
    enriched = _drop_overlap(
        enriched,
        [c for c in div_df.columns if c not in {"compound_id", "smiles", "label"}],
    )
    enriched = enriched.merge(
        div_df.drop(columns=[c for c in ("smiles", "label") if c in div_df.columns]),
        on="compound_id",
        how="left",
    )

    enriched_path = output_dir / "phase2_enriched.csv"
    subset_path = output_dir / "diverse_subset.csv"
    scaffold_path = output_dir / "scaffold_summary.csv"
    audit_path = output_dir / "audit_log.json"
    enriched.to_csv(enriched_path, index=False)
    if "sphere_selected" in enriched.columns:
        selected_mask = enriched["sphere_selected"].astype(str).str.lower().isin({"true", "1"})
        enriched.loc[selected_mask].to_csv(subset_path, index=False)
    else:
        pd.DataFrame().to_csv(subset_path, index=False)
    write_csv(scaffold_path, scaffold_summary_rows(diversity_records))

    summary = {
        "n_input": int(len(frame)),
        "n_reference_actives": novelty_summary.n_reference_actives,
        "properties": property_summary.model_dump(mode="json"),
        "novelty": novelty_summary.model_dump(mode="json"),
        "diversity": diversity_summary.model_dump(mode="json"),
    }
    manifest = build_audit(
        command="phase2",
        config={
            "properties": properties_config.model_dump(mode="json"),
            "novelty": novelty_config.model_dump(mode="json"),
            "diversity": diversity_config.model_dump(mode="json"),
        },
        input_path=input_path,
        output_dir=output_dir,
        summary=summary,
        outputs=[
            str(enriched_path.resolve()),
            str(subset_path.resolve()),
            str(scaffold_path.resolve()),
            str(audit_path.resolve()),
        ],
    )
    write_audit(audit_path, manifest)
    return {
        "property_records": property_records,
        "novelty_records": novelty_records,
        "diversity_records": diversity_records,
        "property_summary": property_summary,
        "novelty_summary": novelty_summary,
        "diversity_summary": diversity_summary,
        "references": references,
        "manifest": manifest,
        "enriched": enriched,
    }
