"""Typer CLI for HitPilot hygiene, properties, novelty, and diversity."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import typer
from rich.console import Console
from rich.table import Table

from hitpilot.data.litpcba import KNOWN_TARGETS, LitPcbaError, fetch_litpcba, normalize_target
from hitpilot.diversity import run_diversity_file
from hitpilot.hygiene import run_hygiene_file
from hitpilot.models import (
    DiversityConfig,
    HygieneConfig,
    HygieneSummary,
    NoveltyConfig,
    PropertiesConfig,
)
from hitpilot.novelty import ChemblFetchError, fetch_chembl_target_actives, run_novelty_file
from hitpilot.phase2 import run_phase2_file
from hitpilot.properties import run_properties_file
from hitpilot.version import __version__

app = typer.Typer(
    name="hitpilot",
    help="HitPilot: auditable hit-triage and medicinal-chemistry copilot.",
    no_args_is_help=True,
)
console = Console()


def _print_summary(summary: HygieneSummary) -> None:
    overview = Table(title="Hygiene summary")
    overview.add_column("Metric")
    overview.add_column("Value", justify="right")
    overview.add_row("Total compounds", str(summary.n_total))
    overview.add_row("Passed", str(summary.n_passed))
    overview.add_row("Dropped", str(summary.n_dropped))
    overview.add_row("Pass rate", f"{summary.pass_rate:.2%}")
    console.print(overview)

    reasons = Table(title="Drop reasons")
    reasons.add_column("Reason")
    reasons.add_column("Count", justify="right")
    if summary.counts_by_drop_reason:
        for reason, count in summary.counts_by_drop_reason.items():
            reasons.add_row(reason, str(count))
    else:
        reasons.add_row("(none)", "0")
    console.print(reasons)

    catalogs = Table(title="FilterCatalog matches (compounds)")
    catalogs.add_column("Catalog")
    catalogs.add_column("Count", justify="right")
    if summary.counts_by_catalog:
        for catalog, count in summary.counts_by_catalog.items():
            catalogs.add_row(catalog, str(count))
    else:
        catalogs.add_row("(none)", "0")
    console.print(catalogs)

    if summary.by_label:
        labels = Table(title="By label")
        labels.add_column("Label")
        labels.add_column("Total", justify="right")
        labels.add_column("Passed", justify="right")
        labels.add_column("Dropped", justify="right")
        labels.add_column("With alerts", justify="right")
        labels.add_column("Pass rate", justify="right")
        for label, breakdown in summary.by_label.items():
            labels.add_row(
                label,
                str(breakdown.n_total),
                str(breakdown.n_passed),
                str(breakdown.n_dropped),
                str(breakdown.n_with_alerts),
                f"{breakdown.pass_rate:.2%}",
            )
        console.print(labels)


def _version_callback(value: bool) -> None:
    if value:
        console.print(f"hitpilot {__version__}")
        raise typer.Exit()


@app.callback()
def main(
    version: bool = typer.Option(
        False,
        "--version",
        help="Print the HitPilot version and exit.",
        is_eager=True,
        callback=_version_callback,
    ),
) -> None:
    """HitPilot command-line interface."""


@app.command("fetch-litpcba")
def fetch_litpcba_command(
    target: str = typer.Argument(..., help="LIT-PCBA target name, e.g. ESR1_ant or PPARG."),
    dest: Path = typer.Option(
        Path("data/litpcba"),
        "--dest",
        "-d",
        help="Directory that will contain data/litpcba/<TARGET>/.",
    ),
    url: Optional[str] = typer.Option(
        None,
        "--url",
        help="Override the official Unistra (or mirror) archive URL.",
    ),
    timeout: int = typer.Option(15, "--timeout", help="Connect timeout in seconds."),
) -> None:
    """Download one public LIT-PCBA target and write a unified CSV."""
    name = normalize_target(target)
    if name not in KNOWN_TARGETS:
        console.print(
            f"[yellow]Warning:[/yellow] {target!r} is not in the built-in LIT-PCBA target list. "
            "Download will still be attempted."
        )
    try:
        result = fetch_litpcba(name, dest=dest, url=url, timeout=timeout)
    except LitPcbaError as exc:
        console.print(f"[red]Failed to fetch {name}:[/red] {exc}")
        raise typer.Exit(code=1) from exc

    table = Table(title=f"LIT-PCBA {result['target']}")
    table.add_column("Field")
    table.add_column("Value")
    for key in (
        "target",
        "protein",
        "description",
        "url",
        "csv_path",
        "n_compounds",
        "n_active",
        "n_inactive",
    ):
        table.add_row(key, str(result.get(key, "")))
    console.print(table)


@app.command("hygiene")
def hygiene_command(
    input_file: Path = typer.Argument(
        ...,
        exists=True,
        readable=True,
        help="CSV, TSV, or SMI input.",
    ),
    output_dir: Path = typer.Option(
        ...,
        "--output-dir",
        "-o",
        help="Directory for passed_compounds.csv, all_compounds.csv, and audit_log.json.",
    ),
    smiles_column: Optional[str] = typer.Option(
        None, "--smiles-column", help="SMILES column name."
    ),
    id_column: Optional[str] = typer.Option(
        None, "--id-column", help="Compound ID column name."
    ),
    label_column: Optional[str] = typer.Option(
        None, "--label-column", help="Activity label column."
    ),
    score_column: Optional[str] = typer.Option(
        None,
        "--score-column",
        help="Numeric score used to pick a keeper when structures collide.",
    ),
    score_ascending: bool = typer.Option(
        False,
        "--score-ascending/--score-descending",
        help="Keep the lowest score (e.g. IC50) instead of the highest.",
    ),
    dedup_key: str = typer.Option(
        "canonical_smiles",
        "--dedup-key",
        help="Deduplicate by canonical_smiles or inchikey.",
    ),
    canonicalize_tautomers: bool = typer.Option(
        False,
        "--canonicalize-tautomers/--no-canonicalize-tautomers",
        help="Canonicalize tautomers with RDKit TautomerEnumerator (off by default).",
    ),
    enable_pains: bool = typer.Option(True, "--enable-pains/--disable-pains"),
    enable_brenk: bool = typer.Option(True, "--enable-brenk/--disable-brenk"),
    enable_nih: bool = typer.Option(True, "--enable-nih/--disable-nih"),
    drop_pains: bool = typer.Option(
        True,
        "--drop-pains/--no-drop-pains",
        help="Drop PAINS matches (flag-only when disabled).",
    ),
    drop_brenk: bool = typer.Option(
        True,
        "--drop-brenk/--no-drop-brenk",
        help="Drop Brenk matches (flag-only when disabled).",
    ),
    drop_nih: bool = typer.Option(
        True,
        "--drop-nih/--no-drop-nih",
        help="Drop NIH matches (flag-only when disabled). Default is to drop.",
    ),
) -> None:
    """Standardize compounds, apply PAINS/Brenk/NIH alerts, and write an audit log."""
    if dedup_key not in {"canonical_smiles", "inchikey"}:
        console.print("[red]--dedup-key must be canonical_smiles or inchikey[/red]")
        raise typer.Exit(code=2)
    config = HygieneConfig(
        canonicalize_tautomers=canonicalize_tautomers,
        dedup_key=dedup_key,  # type: ignore[arg-type]
        score_column=score_column,
        score_ascending=score_ascending,
        enable_pains=enable_pains,
        enable_brenk=enable_brenk,
        enable_nih=enable_nih,
        drop_pains=drop_pains,
        drop_brenk=drop_brenk,
        drop_nih=drop_nih,
    )
    manifest = run_hygiene_file(
        input_file,
        output_dir,
        config,
        smiles_column=smiles_column,
        id_column=id_column,
        label_column=label_column,
        score_column=score_column,
    )
    console.print(f"[green]Wrote[/green] {output_dir / 'passed_compounds.csv'}")
    console.print(f"[green]Wrote[/green] {output_dir / 'all_compounds.csv'}")
    console.print(f"[green]Wrote[/green] {output_dir / 'audit_log.json'}")
    console.print(
        f"HitPilot {manifest.hitpilot_version}  RDKit {manifest.rdkit_version}  "
        f"SHA-256 {manifest.input_sha256[:12]}…  {manifest.timestamp_utc}"
    )
    _print_summary(manifest.summary)


def _print_kv_table(title: str, rows: list[tuple[str, object]]) -> None:
    table = Table(title=title)
    table.add_column("Metric")
    table.add_column("Value", justify="right")
    for key, value in rows:
        if isinstance(value, float):
            rendered = f"{value:.3f}" if abs(value) < 1000 else f"{value:.2f}"
        else:
            rendered = str(value)
        table.add_row(key, rendered)
    console.print(table)


@app.command("properties")
def properties_command(
    input_file: Path = typer.Argument(
        ...,
        exists=True,
        readable=True,
        help="Compound table (typically Phase 1 passed_compounds.csv).",
    ),
    output_dir: Path = typer.Option(
        ...,
        "--output-dir",
        "-o",
        help="Directory for properties.csv and audit_log.json.",
    ),
    smiles_column: Optional[str] = typer.Option(None, "--smiles-column"),
    id_column: Optional[str] = typer.Option(None, "--id-column"),
    label_column: Optional[str] = typer.Option(None, "--label-column"),
) -> None:
    """Compute RDKit descriptors, Ertl SA score, and Lipinski/Veber flags."""
    records, summary, manifest = run_properties_file(
        input_file,
        output_dir,
        PropertiesConfig(),
        smiles_column=smiles_column,
        id_column=id_column,
        label_column=label_column,
    )
    del records
    console.print(f"[green]Wrote[/green] {output_dir / 'properties.csv'}")
    console.print(f"[green]Wrote[/green] {output_dir / 'audit_log.json'}")
    console.print(
        f"HitPilot {manifest.hitpilot_version}  RDKit {manifest.rdkit_version}  "
        f"SHA-256 {manifest.input_sha256[:12]}…  {manifest.timestamp_utc}"
    )
    _print_kv_table(
        "Property summary",
        [
            ("Total compounds", summary.n_total),
            ("Valid", summary.n_valid),
            ("Median MW", summary.median_mw),
            ("Median cLogP", summary.median_clogp),
            ("Median TPSA", summary.median_tpsa),
            ("Median QED", summary.median_qed),
            ("Median SA score", summary.median_sa_score),
            ("Lipinski pass", f"{summary.lipinski_pass_count} ({summary.lipinski_pass_rate:.2%})"),
            (
                "Lipinski strict pass",
                f"{summary.lipinski_strict_pass_count} ({summary.lipinski_strict_pass_rate:.2%})",
            ),
            ("Veber pass", f"{summary.veber_pass_count} ({summary.veber_pass_rate:.2%})"),
        ],
    )
    if summary.by_label:
        table = Table(title="By label")
        table.add_column("Label")
        table.add_column("N", justify="right")
        table.add_column("Med MW", justify="right")
        table.add_column("Med cLogP", justify="right")
        table.add_column("Med QED", justify="right")
        table.add_column("Lipinski", justify="right")
        table.add_column("Veber", justify="right")
        for label, stats in summary.by_label.items():
            table.add_row(
                label,
                str(stats.n_valid),
                f"{stats.median_mw:.1f}" if stats.median_mw is not None else "",
                f"{stats.median_clogp:.2f}" if stats.median_clogp is not None else "",
                f"{stats.median_qed:.3f}" if stats.median_qed is not None else "",
                f"{stats.lipinski_pass_rate:.2%}",
                f"{stats.veber_pass_rate:.2%}",
            )
        console.print(table)


@app.command("fetch-chembl")
def fetch_chembl_command(
    target_chembl_id: str = typer.Argument("CHEMBL206", help="ChEMBL target ID."),
    pchembl_min: float = typer.Option(6.0, "--pchembl-min", help="Minimum pChEMBL value."),
    max_records: int = typer.Option(2000, "--max-records", help="Max activity records to pull."),
    cache_dir: Path = typer.Option(Path("data/chembl"), "--cache-dir", help="Cache directory."),
    refresh: bool = typer.Option(False, "--refresh", help="Ignore any existing cache."),
    timeout: int = typer.Option(15, "--timeout", help="Connect timeout in seconds."),
) -> None:
    """Download unique ChEMBL actives for a target and cache them as CSV."""
    try:
        actives = fetch_chembl_target_actives(
            target_chembl_id=target_chembl_id,
            pchembl_min=pchembl_min,
            max_records=max_records,
            cache_dir=cache_dir,
            refresh_cache=refresh,
            timeout=timeout,
        )
    except ChemblFetchError as exc:
        console.print(f"[red]Failed to fetch {target_chembl_id}:[/red] {exc}")
        raise typer.Exit(code=1) from exc
    csv_path = cache_dir / f"{target_chembl_id.upper()}_actives.csv"
    pchembls = [item.pchembl_value for item in actives]
    _print_kv_table(
        f"ChEMBL {target_chembl_id} actives",
        [
            ("Unique actives", len(actives)),
            ("pChEMBL threshold", pchembl_min),
            ("Max activity records", max_records),
            ("Min pChEMBL kept", min(pchembls) if pchembls else ""),
            ("Max pChEMBL kept", max(pchembls) if pchembls else ""),
            ("Cache", str(csv_path)),
        ],
    )


@app.command("novelty")
def novelty_command(
    input_file: Path = typer.Argument(..., exists=True, readable=True),
    output_dir: Path = typer.Option(..., "--output-dir", "-o"),
    target_chembl_id: str = typer.Option("CHEMBL206", "--target-chembl-id"),
    pchembl_min: float = typer.Option(6.0, "--pchembl-min"),
    max_records: int = typer.Option(2000, "--max-records"),
    cache_dir: Path = typer.Option(Path("data/chembl"), "--cache-dir"),
    refresh: bool = typer.Option(False, "--refresh"),
    smiles_column: Optional[str] = typer.Option(None, "--smiles-column"),
    id_column: Optional[str] = typer.Option(None, "--id-column"),
    label_column: Optional[str] = typer.Option(None, "--label-column"),
) -> None:
    """Score Morgan (ECFP4-like) novelty against ChEMBL reference actives."""
    config = NoveltyConfig(
        target_chembl_id=target_chembl_id,
        pchembl_min=pchembl_min,
        max_records=max_records,
        cache_dir=str(cache_dir),
        refresh_cache=refresh,
    )
    try:
        records, summary, manifest = run_novelty_file(
            input_file,
            output_dir,
            config,
            smiles_column=smiles_column,
            id_column=id_column,
            label_column=label_column,
        )
    except ChemblFetchError as exc:
        console.print(f"[red]Failed to load ChEMBL references:[/red] {exc}")
        raise typer.Exit(code=1) from exc
    del records
    console.print(f"[green]Wrote[/green] {output_dir / 'novelty.csv'}")
    console.print(f"[green]Wrote[/green] {output_dir / 'audit_log.json'}")
    console.print(
        f"HitPilot {manifest.hitpilot_version}  RDKit {manifest.rdkit_version}  "
        f"SHA-256 {manifest.input_sha256[:12]}…  {manifest.timestamp_utc}"
    )
    rows = [
        ("Total compounds", summary.n_total),
        ("Reference actives", summary.n_reference_actives),
        ("Median max Tanimoto", summary.median_max_tanimoto),
    ]
    rows.extend((tier, count) for tier, count in summary.counts_by_tier.items())
    _print_kv_table("Novelty summary", rows)


@app.command("diversity")
def diversity_command(
    input_file: Path = typer.Argument(..., exists=True, readable=True),
    output_dir: Path = typer.Option(..., "--output-dir", "-o"),
    tanimoto_cutoff: float = typer.Option(
        0.60, "--tanimoto-cutoff", help="Butina similarity cutoff."
    ),
    exclusion_similarity: float = typer.Option(
        0.65, "--exclusion-similarity", help="Sphere-exclusion Tanimoto threshold."
    ),
    score_column: Optional[str] = typer.Option(None, "--score-column"),
    score_ascending: bool = typer.Option(False, "--score-ascending/--score-descending"),
    rank_by_score: bool = typer.Option(False, "--rank-by-score/--input-order"),
    smiles_column: Optional[str] = typer.Option(None, "--smiles-column"),
    id_column: Optional[str] = typer.Option(None, "--id-column"),
    label_column: Optional[str] = typer.Option(None, "--label-column"),
) -> None:
    """Assign Bemis-Murcko scaffolds, Butina clusters, and a diverse subset."""
    config = DiversityConfig(
        tanimoto_cutoff=tanimoto_cutoff,
        exclusion_similarity=exclusion_similarity,
        score_column=score_column,
        score_ascending=score_ascending,
        rank_by_score=rank_by_score or bool(score_column),
    )
    records, summary, manifest = run_diversity_file(
        input_file,
        output_dir,
        config,
        smiles_column=smiles_column,
        id_column=id_column,
        label_column=label_column,
        score_column=score_column,
    )
    del records
    console.print(f"[green]Wrote[/green] {output_dir / 'diversity.csv'}")
    console.print(f"[green]Wrote[/green] {output_dir / 'diverse_subset.csv'}")
    console.print(f"[green]Wrote[/green] {output_dir / 'scaffold_summary.csv'}")
    console.print(
        f"HitPilot {manifest.hitpilot_version}  RDKit {manifest.rdkit_version}  "
        f"SHA-256 {manifest.input_sha256[:12]}…  {manifest.timestamp_utc}"
    )
    _print_kv_table(
        "Diversity summary",
        [
            ("Total compounds", summary.n_total),
            ("Unique Murcko scaffolds", summary.n_unique_murcko_scaffolds),
            ("Unique generic scaffolds", summary.n_unique_generic_scaffolds),
            ("Acyclic", summary.n_acyclic),
            ("Butina clusters", summary.n_clusters),
            ("Top cluster sizes", ", ".join(str(s) for s in summary.top_cluster_sizes[:5])),
            ("Sphere-selected", summary.n_sphere_selected),
            (
                "Actives selected",
                (
                    f"{summary.n_actives_selected}/{summary.n_actives_total}"
                    if summary.n_actives_total is not None
                    else ""
                ),
            ),
        ],
    )


@app.command("phase2")
def phase2_command(
    input_file: Path = typer.Argument(
        ...,
        exists=True,
        readable=True,
        help="Phase 1 passed_compounds.csv or another compound table.",
    ),
    output_dir: Path = typer.Option(
        ...,
        "--output-dir",
        "-o",
        help="Directory for phase2_enriched.csv and related outputs.",
    ),
    target_chembl_id: str = typer.Option("CHEMBL206", "--target-chembl-id"),
    pchembl_min: float = typer.Option(6.0, "--pchembl-min"),
    max_records: int = typer.Option(2000, "--max-records"),
    cache_dir: Path = typer.Option(Path("data/chembl"), "--cache-dir"),
    refresh: bool = typer.Option(False, "--refresh"),
    tanimoto_cutoff: float = typer.Option(0.60, "--tanimoto-cutoff"),
    exclusion_similarity: float = typer.Option(0.65, "--exclusion-similarity"),
    score_column: Optional[str] = typer.Option(None, "--score-column"),
    smiles_column: Optional[str] = typer.Option(None, "--smiles-column"),
    id_column: Optional[str] = typer.Option(None, "--id-column"),
    label_column: Optional[str] = typer.Option(None, "--label-column"),
) -> None:
    """Run properties, novelty, and diversity and write a unified Phase 2 dossier."""
    novelty_config = NoveltyConfig(
        target_chembl_id=target_chembl_id,
        pchembl_min=pchembl_min,
        max_records=max_records,
        cache_dir=str(cache_dir),
        refresh_cache=refresh,
    )
    diversity_config = DiversityConfig(
        tanimoto_cutoff=tanimoto_cutoff,
        exclusion_similarity=exclusion_similarity,
        score_column=score_column,
        rank_by_score=bool(score_column),
    )
    try:
        result = run_phase2_file(
            input_file,
            output_dir,
            properties_config=PropertiesConfig(),
            novelty_config=novelty_config,
            diversity_config=diversity_config,
            smiles_column=smiles_column,
            id_column=id_column,
            label_column=label_column,
            score_column=score_column,
        )
    except ChemblFetchError as exc:
        console.print(f"[red]Failed to load ChEMBL references:[/red] {exc}")
        raise typer.Exit(code=1) from exc
    manifest = result["manifest"]
    props = result["property_summary"]
    nov = result["novelty_summary"]
    div = result["diversity_summary"]
    console.print(f"[green]Wrote[/green] {output_dir / 'phase2_enriched.csv'}")
    console.print(f"[green]Wrote[/green] {output_dir / 'diverse_subset.csv'}")
    console.print(f"[green]Wrote[/green] {output_dir / 'scaffold_summary.csv'}")
    console.print(f"[green]Wrote[/green] {output_dir / 'audit_log.json'}")
    console.print(
        f"HitPilot {manifest.hitpilot_version}  RDKit {manifest.rdkit_version}  "
        f"SHA-256 {manifest.input_sha256[:12]}…  {manifest.timestamp_utc}"
    )
    _print_kv_table(
        "Phase 2 properties",
        [
            ("Compounds", props.n_total),
            ("Median MW", props.median_mw),
            ("Median cLogP", props.median_clogp),
            ("Median TPSA", props.median_tpsa),
            ("Median QED", props.median_qed),
            ("Median SA score", props.median_sa_score),
            ("Lipinski pass", f"{props.lipinski_pass_count} ({props.lipinski_pass_rate:.2%})"),
            ("Veber pass", f"{props.veber_pass_count} ({props.veber_pass_rate:.2%})"),
        ],
    )
    nov_rows = [
        ("Reference actives", nov.n_reference_actives),
        ("Median max Tanimoto", nov.median_max_tanimoto),
    ]
    nov_rows.extend((tier, count) for tier, count in nov.counts_by_tier.items())
    _print_kv_table("Phase 2 novelty", nov_rows)
    _print_kv_table(
        "Phase 2 diversity",
        [
            ("Unique Murcko scaffolds", div.n_unique_murcko_scaffolds),
            ("Unique generic scaffolds", div.n_unique_generic_scaffolds),
            ("Butina clusters", div.n_clusters),
            ("Sphere-selected", div.n_sphere_selected),
            (
                "Actives captured",
                (
                    f"{div.n_actives_selected}/{div.n_actives_total}"
                    if div.n_actives_total is not None
                    else ""
                ),
            ),
        ],
    )


if __name__ == "__main__":
    app()
