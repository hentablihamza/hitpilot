"""Typer CLI for HitPilot Phase 1."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import typer
from rich.console import Console
from rich.table import Table

from hitpilot.data.litpcba import KNOWN_TARGETS, LitPcbaError, fetch_litpcba, normalize_target
from hitpilot.hygiene import run_hygiene_file
from hitpilot.models import HygieneConfig, HygieneSummary
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


if __name__ == "__main__":
    app()
