"""Shared table I/O, hashing, and audit-log helpers."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import pandas as pd

from hitpilot.version import __version__

_SMILES_CANDIDATES = (
    "standardized_smiles",
    "canonical_smiles",
    "smiles",
    "SMILES",
    "Smiles",
    "input_smiles",
    "iso_smiles",
    "mol",
)
_ID_CANDIDATES = (
    "compound_id",
    "compound",
    "id",
    "ID",
    "name",
    "Name",
    "molecule_chembl_id",
    "chembl_id",
    "zinc_id",
)
_LABEL_CANDIDATES = ("label", "Label", "activity", "active", "class", "y")


def utc_timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def rdkit_version() -> str:
    return getattr(__import__("rdkit"), "__version__", "unknown")


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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


def load_compound_frame(
    path: str | Path,
    *,
    smiles_column: str | None = None,
    id_column: str | None = None,
    label_column: str | None = None,
    score_column: str | None = None,
) -> pd.DataFrame:
    """Load a compound table and expose a resolved ``smiles`` column.

    Prefers already-standardized SMILES when present (Phase 1 output),
    otherwise falls back to ``smiles`` / ``input_smiles``. All original
    columns are preserved.
    """
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix in {".tsv", ".txt"}:
        frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    elif suffix in {".smi", ".smiles"}:
        from hitpilot.hygiene import read_compound_table

        frame = read_compound_table(path)
        frame = frame.astype(str)
    else:
        frame = pd.read_csv(path, dtype=str, keep_default_na=False)

    smiles_col = _pick_column(frame.columns, _SMILES_CANDIDATES, smiles_column)
    if smiles_col is None:
        raise ValueError(f"No SMILES column found in {path}. Columns: {list(frame.columns)}")
    id_col = _pick_column(frame.columns, _ID_CANDIDATES, id_column)
    label_col = _pick_column(frame.columns, _LABEL_CANDIDATES, label_column)
    score_col = _pick_column(
        frame.columns,
        ("score", "Score", "docking_score", "affinity"),
        score_column,
    )

    out = frame.copy()
    if "compound_id" not in out.columns:
        if id_col is not None:
            out["compound_id"] = out[id_col].astype(str)
        else:
            out["compound_id"] = [f"cmpd_{i + 1:06d}" for i in range(len(out))]
    else:
        out["compound_id"] = out["compound_id"].astype(str)
    # Per-row: use the resolved smiles column, but if a row is empty fall back.
    resolved: list[str] = []
    for _, row in out.iterrows():
        raw = row[smiles_col] if smiles_col in row else ""
        value = str(raw).strip()
        if not value:
            for alt in _SMILES_CANDIDATES:
                if alt in out.columns and str(row[alt]).strip():
                    value = str(row[alt]).strip()
                    break
        resolved.append(value)
    out["smiles"] = resolved
    if label_col is not None and "label" not in out.columns:
        out["label"] = out[label_col]
    if score_col is not None and "score" not in out.columns:
        out["score"] = out[score_col]
    return out


def write_csv(
    path: str | Path,
    rows: list[dict[str, Any]],
    fieldnames: list[str] | None = None,
) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        empty = pd.DataFrame(columns=fieldnames or [])
        empty.to_csv(path, index=False)
        return path
    frame = pd.DataFrame(rows)
    if fieldnames:
        for name in fieldnames:
            if name not in frame.columns:
                frame[name] = ""
        frame = frame[fieldnames]
    frame.to_csv(path, index=False)
    return path


def write_audit(path: str | Path, payload: Any) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if hasattr(payload, "model_dump"):
        data = payload.model_dump(mode="json")
    else:
        data = payload
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    return path


def build_audit(
    *,
    command: str,
    config: Any,
    input_path: str | Path,
    output_dir: str | Path,
    summary: Any,
    outputs: list[str] | None = None,
):
    from hitpilot.models import PhaseAuditManifest

    config_data = config.model_dump(mode="json") if hasattr(config, "model_dump") else config
    summary_data = summary.model_dump(mode="json") if hasattr(summary, "model_dump") else summary
    return PhaseAuditManifest(
        timestamp_utc=utc_timestamp(),
        hitpilot_version=__version__,
        rdkit_version=rdkit_version(),
        command=command,
        config=config_data,
        input_path=str(Path(input_path).resolve()),
        input_sha256=sha256_file(input_path),
        output_dir=str(Path(output_dir).resolve()),
        summary=summary_data,
        outputs=outputs or [],
    )
