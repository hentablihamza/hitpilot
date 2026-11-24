"""RDKit structure hygiene, FilterCatalog alerts, and auditable batch I/O."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import pandas as pd
from rdkit import Chem, rdBase
from rdkit.Chem import FilterCatalog
from rdkit.Chem.MolStandardize import rdMolStandardize

from hitpilot.models import (
    AlertHit,
    AuditManifest,
    HygieneConfig,
    HygieneRecord,
    HygieneSummary,
    LabelBreakdown,
)
from hitpilot.version import __version__

DROP_INVALID = "invalid_smiles"
DROP_SANITIZE = "sanitization_failed"
DROP_EMPTY = "empty_after_strip"
DROP_INORGANIC = "no_organic_carbon"
DROP_DUPLICATE = "duplicate_structure"
DROP_PAINS = "pains_alert"
DROP_BRENK = "brenk_alert"
DROP_NIH = "nih_alert"

STRUCTURE_FAILURES = {DROP_INVALID, DROP_SANITIZE, DROP_EMPTY, DROP_INORGANIC}

_FILTERSET_TO_CATALOG = {
    "PAINS_A": "PAINS_A",
    "PAINS_B": "PAINS_B",
    "PAINS_C": "PAINS_C",
    "PAINS": "PAINS",
    "Brenk": "BRENK",
    "BRENK": "BRENK",
    "NIH": "NIH",
}

_SMILES_COLUMNS = ("smiles", "SMILES", "Smiles", "canonical_smiles", "iso_smiles", "mol")
_ID_COLUMNS = (
    "compound_id",
    "compound",
    "id",
    "ID",
    "name",
    "Name",
    "zinc_id",
    "chembl_id",
    "ligand",
)
_LABEL_COLUMNS = ("label", "Label", "activity", "active", "class", "y")


@dataclass
class StandardizeResult:
    mol: Chem.Mol | None = None
    smiles: str | None = None
    inchikey: str | None = None
    drop_reasons: list[str] = field(default_factory=list)
    error: str | None = None


def _silence_rdkit_logs() -> None:
    rdBase.DisableLog("rdApp.error")
    rdBase.DisableLog("rdApp.warning")


def build_filter_catalog(config: HygieneConfig) -> FilterCatalog.FilterCatalog | None:
    """Assemble an RDKit FilterCatalog from the enabled PAINS / Brenk / NIH sets."""
    catalogs = FilterCatalog.FilterCatalogParams.FilterCatalogs
    params = FilterCatalog.FilterCatalogParams()
    added = False
    if config.enable_pains:
        params.AddCatalog(catalogs.PAINS_A)
        params.AddCatalog(catalogs.PAINS_B)
        params.AddCatalog(catalogs.PAINS_C)
        added = True
    if config.enable_brenk:
        params.AddCatalog(catalogs.BRENK)
        added = True
    if config.enable_nih:
        params.AddCatalog(catalogs.NIH)
        added = True
    if not added:
        return None
    return FilterCatalog.FilterCatalog(params)


def _has_organic_carbon(mol: Chem.Mol) -> bool:
    return any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())


def standardize_smiles(smiles: str, config: HygieneConfig | None = None) -> StandardizeResult:
    """Parse, sanitize, strip salts, neutralize, and canonicalize one SMILES string."""
    config = config or HygieneConfig()
    raw = "" if smiles is None else str(smiles).strip()
    if not raw:
        return StandardizeResult(drop_reasons=[DROP_INVALID], error="empty SMILES")

    mol = Chem.MolFromSmiles(raw, sanitize=False)
    if mol is None:
        return StandardizeResult(drop_reasons=[DROP_INVALID], error="MolFromSmiles returned None")

    try:
        flag = Chem.SanitizeMol(mol, catchErrors=True)
    except Exception as exc:  # noqa: BLE001 - RDKit can raise several valence errors
        return StandardizeResult(
            drop_reasons=[DROP_SANITIZE],
            error=f"{type(exc).__name__}: {exc}",
        )
    if flag != Chem.SanitizeFlags.SANITIZE_NONE:
        return StandardizeResult(
            drop_reasons=[DROP_SANITIZE],
            error=f"SanitizeMol failed with flags={int(flag)}",
        )

    try:
        chooser = rdMolStandardize.LargestFragmentChooser(preferOrganic=True)
        mol = chooser.choose(mol)
        if mol is None or mol.GetNumAtoms() == 0:
            return StandardizeResult(
                drop_reasons=[DROP_EMPTY],
                error="no atoms remaining after fragment selection",
            )
        if not _has_organic_carbon(mol):
            return StandardizeResult(
                mol=mol,
                drop_reasons=[DROP_INORGANIC],
                error="no organic carbon after fragment selection",
            )
        normalizer = rdMolStandardize.Normalizer()
        mol = normalizer.normalize(mol)
        uncharger = rdMolStandardize.Uncharger()
        mol = uncharger.uncharge(mol)
        if config.canonicalize_tautomers:
            enumerator = rdMolStandardize.TautomerEnumerator()
            mol = enumerator.Canonicalize(mol)
        Chem.SanitizeMol(mol)
    except Exception as exc:  # noqa: BLE001
        return StandardizeResult(
            drop_reasons=[DROP_SANITIZE],
            error=f"{type(exc).__name__}: {exc}",
        )

    if mol is None or mol.GetNumAtoms() == 0:
        return StandardizeResult(
            drop_reasons=[DROP_EMPTY],
            error="empty structure after standardization",
        )
    if not _has_organic_carbon(mol):
        return StandardizeResult(
            mol=mol,
            drop_reasons=[DROP_INORGANIC],
            error="no organic carbon after standardization",
        )

    canon = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
    inchikey: str | None = None
    try:
        inchikey = Chem.MolToInchiKey(mol)
    except Exception as exc:  # noqa: BLE001
        inchikey = None
        # Keep the standardized structure; InChIKey failure is not a drop reason.
        _ = exc
    return StandardizeResult(mol=mol, smiles=canon, inchikey=inchikey)


def collect_alerts(mol: Chem.Mol, catalog: FilterCatalog.FilterCatalog | None) -> list[AlertHit]:
    """Return every FilterCatalog match on a standardized molecule."""
    if catalog is None or mol is None:
        return []
    hits: list[AlertHit] = []
    seen: set[tuple[str, str]] = set()
    for entry in catalog.GetMatches(mol):
        props = set(entry.GetPropList())
        raw_set = entry.GetProp("FilterSet") if "FilterSet" in props else ""
        catalog_name = _FILTERSET_TO_CATALOG.get(raw_set, raw_set.upper() or "UNKNOWN")
        if entry.GetDescription():
            rule = entry.GetDescription()
        elif "description" in props:
            rule = entry.GetProp("description")
        else:
            rule = ""
        key = (catalog_name, rule)
        if key in seen:
            continue
        seen.add(key)
        hits.append(AlertHit(catalog=catalog_name, rule=rule))
    return hits


def _alert_drop_reasons(alerts: Iterable[AlertHit], config: HygieneConfig) -> list[str]:
    reasons: list[str] = []
    families = {hit.catalog for hit in alerts}
    has_pains = any(name.startswith("PAINS") for name in families)
    if has_pains and config.enable_pains and config.drop_pains:
        reasons.append(DROP_PAINS)
    if "BRENK" in families and config.enable_brenk and config.drop_brenk:
        reasons.append(DROP_BRENK)
    if "NIH" in families and config.enable_nih and config.drop_nih:
        reasons.append(DROP_NIH)
    return reasons


def _parse_score(value: Any) -> float | None:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    text = str(value).strip()
    if text == "" or text.lower() in {"nan", "none", "na"}:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def process_compound(
    compound_id: str,
    smiles: str,
    config: HygieneConfig,
    catalog: FilterCatalog.FilterCatalog | None,
    *,
    label: str | None = None,
    score: float | None = None,
    input_index: int = 0,
) -> HygieneRecord:
    """Standardize one compound and attach FilterCatalog alerts."""
    result = standardize_smiles(smiles, config)
    record = HygieneRecord(
        compound_id=str(compound_id),
        input_smiles="" if smiles is None else str(smiles),
        standardized_smiles=result.smiles,
        inchikey=result.inchikey,
        drop_reasons=list(result.drop_reasons),
        label=label,
        score=score,
        error=result.error,
        input_index=input_index,
    )
    if result.drop_reasons or result.mol is None:
        record.passed = False
        return record
    record.alerts = collect_alerts(result.mol, catalog)
    record.drop_reasons.extend(_alert_drop_reasons(record.alerts, config))
    record.passed = not record.drop_reasons
    return record


def _dedup_key(record: HygieneRecord, config: HygieneConfig) -> str | None:
    if config.dedup_key == "inchikey":
        if record.inchikey:
            return f"inchikey:{record.inchikey}"
        if record.standardized_smiles:
            return f"smiles:{record.standardized_smiles}"
        return None
    if record.standardized_smiles:
        return f"smiles:{record.standardized_smiles}"
    return None


def _is_better(candidate: HygieneRecord, current: HygieneRecord, config: HygieneConfig) -> bool:
    """Return True if candidate should replace the current keeper."""
    cand_score, curr_score = candidate.score, current.score
    if cand_score is not None and curr_score is not None and cand_score != curr_score:
        if config.score_ascending:
            return cand_score < curr_score
        return cand_score > curr_score
    if cand_score is not None and curr_score is None:
        return True
    if cand_score is None and curr_score is not None:
        return False
    return candidate.input_index < current.input_index


def deduplicate_records(records: list[HygieneRecord], config: HygieneConfig) -> list[HygieneRecord]:
    """Keep one standardized structure per key; tag the rest as duplicates."""
    groups: dict[str, list[int]] = defaultdict(list)
    for idx, record in enumerate(records):
        if set(record.drop_reasons) & STRUCTURE_FAILURES:
            continue
        key = _dedup_key(record, config)
        if key is None:
            continue
        groups[key].append(idx)

    for indices in groups.values():
        if len(indices) < 2:
            continue
        keeper_idx = indices[0]
        for idx in indices[1:]:
            if _is_better(records[idx], records[keeper_idx], config):
                keeper_idx = idx
        keeper_id = records[keeper_idx].compound_id
        for idx in indices:
            if idx == keeper_idx:
                continue
            record = records[idx]
            if DROP_DUPLICATE not in record.drop_reasons:
                # Keep structure-quality reasons first, then duplicate, then alerts.
                structure = [r for r in record.drop_reasons if r in STRUCTURE_FAILURES]
                alerts = [r for r in record.drop_reasons if r not in STRUCTURE_FAILURES]
                record.drop_reasons = structure + [DROP_DUPLICATE] + alerts
            record.duplicate_of = keeper_id
            record.passed = False
    return records


def summarize_records(records: list[HygieneRecord]) -> HygieneSummary:
    n_total = len(records)
    n_passed = sum(1 for rec in records if rec.passed)
    reason_counts: dict[str, int] = defaultdict(int)
    catalog_counts: dict[str, int] = defaultdict(int)
    for rec in records:
        for reason in rec.drop_reasons:
            reason_counts[reason] += 1
        for catalog in {hit.catalog for hit in rec.alerts}:
            catalog_counts[catalog] += 1

    by_label: dict[str, LabelBreakdown] = {}
    labels = sorted({rec.label for rec in records if rec.label})
    for label in labels:
        subset = [rec for rec in records if rec.label == label]
        label_reasons: dict[str, int] = defaultdict(int)
        for rec in subset:
            for reason in rec.drop_reasons:
                label_reasons[reason] += 1
        n_sub = len(subset)
        n_pass = sum(1 for rec in subset if rec.passed)
        by_label[label] = LabelBreakdown(
            n_total=n_sub,
            n_passed=n_pass,
            n_dropped=n_sub - n_pass,
            n_with_alerts=sum(1 for rec in subset if rec.alerts),
            pass_rate=(n_pass / n_sub) if n_sub else 0.0,
            counts_by_drop_reason=dict(sorted(label_reasons.items())),
        )

    return HygieneSummary(
        n_total=n_total,
        n_passed=n_passed,
        n_dropped=n_total - n_passed,
        pass_rate=(n_passed / n_total) if n_total else 0.0,
        counts_by_drop_reason=dict(sorted(reason_counts.items())),
        counts_by_catalog=dict(sorted(catalog_counts.items())),
        by_label=by_label,
    )


def run_hygiene(
    rows: Iterable[dict[str, Any]],
    config: HygieneConfig | None = None,
    *,
    smiles_field: str = "smiles",
    id_field: str = "compound_id",
    label_field: str = "label",
    score_field: str | None = None,
) -> tuple[list[HygieneRecord], HygieneSummary]:
    """Run standardization, alerts, and deduplication on in-memory rows."""
    _silence_rdkit_logs()
    config = config or HygieneConfig()
    catalog = build_filter_catalog(config)
    score_name = score_field or config.score_column
    records: list[HygieneRecord] = []
    for index, row in enumerate(rows):
        compound_id = row.get(id_field) or row.get("compound_id") or f"cmpd_{index + 1:06d}"
        smiles = row.get(smiles_field, row.get("smiles", ""))
        raw_label = row.get(label_field, row.get("label"))
        if raw_label is None or (isinstance(raw_label, float) and pd.isna(raw_label)):
            label = None
        else:
            label = str(raw_label)
        if label == "":
            label = None
        score = _parse_score(row.get(score_name)) if score_name else _parse_score(row.get("score"))
        records.append(
            process_compound(
                str(compound_id),
                "" if smiles is None else str(smiles),
                config,
                catalog,
                label=label,
                score=score,
                input_index=index,
            )
        )
    records = deduplicate_records(records, config)
    # Recompute passed now that duplicates are tagged (alerts already applied).
    for rec in records:
        rec.passed = not rec.drop_reasons
    return records, summarize_records(records)


def _pick_column(
    columns: Iterable[str],
    candidates: tuple[str, ...],
    explicit: str | None,
) -> str | None:
    if explicit:
        if explicit in columns:
            return explicit
        lowered = {col.lower(): col for col in columns}
        return lowered.get(explicit.lower(), explicit)
    for name in candidates:
        if name in columns:
            return name
    lowered = {col.lower(): col for col in columns}
    for name in candidates:
        if name.lower() in lowered:
            return lowered[name.lower()]
    return None


def _read_smi(path: Path) -> pd.DataFrame:
    ids: list[str] = []
    smiles: list[str] = []
    extras: list[str] = []
    with path.open(newline="", encoding="utf-8") as handle:
        for index, line in enumerate(handle):
            text = line.strip()
            if not text or text.startswith("#"):
                continue
            parts = text.replace(",", " ").replace("\t", " ").split()
            if not parts:
                continue
            # Header-like first line
            if index == 0 and parts[0].lower() in {"smiles", "canonical_smiles"}:
                continue
            smiles.append(parts[0])
            if len(parts) >= 2:
                ids.append(parts[1])
            else:
                ids.append(f"cmpd_{len(ids) + 1:06d}")
            extras.append(" ".join(parts[2:]) if len(parts) > 2 else "")
    frame = pd.DataFrame({"compound_id": ids, "smiles": smiles})
    if any(extras):
        frame["extra"] = extras
    return frame


def read_compound_table(
    path: str | Path,
    *,
    smiles_column: str | None = None,
    id_column: str | None = None,
    label_column: str | None = None,
    score_column: str | None = None,
) -> pd.DataFrame:
    """Load a CSV, TSV, or SMI file into a table with canonical column names."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    suffix = path.suffix.lower()
    if suffix in {".smi", ".smiles"}:
        frame = _read_smi(path)
    elif suffix in {".tsv", ".txt"}:
        frame = pd.read_csv(path, sep="\t")
    elif suffix == ".csv":
        frame = pd.read_csv(path)
    else:
        # Peek at the first line to guess the delimiter.
        sample = path.read_text(encoding="utf-8", errors="replace").splitlines()
        first = next((line for line in sample if line.strip() and not line.startswith("#")), "")
        if "\t" in first:
            frame = pd.read_csv(path, sep="\t")
        elif "," in first:
            frame = pd.read_csv(path)
        else:
            frame = _read_smi(path)

    smiles_col = _pick_column(frame.columns, _SMILES_COLUMNS, smiles_column)
    if smiles_col is None:
        raise ValueError(f"No SMILES column found in {path}. Columns: {list(frame.columns)}")
    id_col = _pick_column(frame.columns, _ID_COLUMNS, id_column)
    label_col = _pick_column(frame.columns, _LABEL_COLUMNS, label_column)
    score_col = _pick_column(
        frame.columns,
        ("score", "Score", "docking_score", "affinity"),
        score_column,
    )

    out = pd.DataFrame(
        {
            "compound_id": (
                frame[id_col].astype(str)
                if id_col is not None
                else [f"cmpd_{i + 1:06d}" for i in range(len(frame))]
            ),
            "smiles": frame[smiles_col].astype(str),
        }
    )
    if label_col is not None:
        out["label"] = frame[label_col]
    if score_col is not None:
        out["score"] = frame[score_col]
    return out


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_records_csv(path: Path, records: list[HygieneRecord]) -> None:
    fieldnames = [
        "compound_id",
        "input_smiles",
        "standardized_smiles",
        "inchikey",
        "passed",
        "drop_reasons",
        "alerts",
        "duplicate_of",
        "label",
        "score",
        "error",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for record in records:
            writer.writerow(record.to_csv_row())


def run_hygiene_file(
    input_path: str | Path,
    output_dir: str | Path,
    config: HygieneConfig | None = None,
    *,
    smiles_column: str | None = None,
    id_column: str | None = None,
    label_column: str | None = None,
    score_column: str | None = None,
) -> AuditManifest:
    """Run hygiene on a CSV/TSV/SMI file and write CSV plus audit-log outputs."""
    input_path = Path(input_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    config = config or HygieneConfig()
    if score_column and not config.score_column:
        config = config.model_copy(update={"score_column": score_column})

    frame = read_compound_table(
        input_path,
        smiles_column=smiles_column,
        id_column=id_column,
        label_column=label_column,
        score_column=score_column or config.score_column,
    )
    records, summary = run_hygiene(
        frame.to_dict(orient="records"),
        config,
        smiles_field="smiles",
        id_field="compound_id",
        label_field="label",
        score_field="score" if "score" in frame.columns else None,
    )

    all_path = output_dir / "all_compounds.csv"
    passed_path = output_dir / "passed_compounds.csv"
    audit_path = output_dir / "audit_log.json"
    _write_records_csv(all_path, records)
    _write_records_csv(passed_path, [rec for rec in records if rec.passed])

    manifest = AuditManifest(
        timestamp_utc=datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        hitpilot_version=__version__,
        rdkit_version=getattr(__import__("rdkit"), "__version__", "unknown"),
        config=config,
        input_path=str(input_path.resolve()),
        input_sha256=_sha256_file(input_path),
        output_dir=str(output_dir.resolve()),
        summary=summary,
    )
    audit_path.write_text(
        json.dumps(manifest.model_dump(mode="json"), indent=2) + "\n",
        encoding="utf-8",
    )
    return manifest
