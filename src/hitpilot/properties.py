"""RDKit physicochemical descriptors, SA score, and Lipinski/Veber flags."""

from __future__ import annotations

import importlib.util
import sys
from collections import defaultdict
from pathlib import Path
from statistics import median
from typing import Any, Iterable

from rdkit import Chem, RDConfig
from rdkit.Chem import QED, Crippen, Descriptors, rdMolDescriptors

from hitpilot.io_utils import build_audit, load_compound_frame, write_audit, write_csv
from hitpilot.models import (
    PropertiesConfig,
    PropertiesSummary,
    PropertyLabelStats,
    PropertyRecord,
)

_SASCORER = None


def load_sascorer():
    """Load RDKit's built-in Ertl SA_Score implementation (never approximate)."""
    global _SASCORER
    if _SASCORER is not None:
        return _SASCORER
    path = Path(RDConfig.RDContribDir) / "SA_Score" / "sascorer.py"
    if not path.is_file():
        raise FileNotFoundError(f"RDKit SA_Score module not found at {path}")
    if "sascorer" in sys.modules:
        _SASCORER = sys.modules["sascorer"]
        return _SASCORER
    spec = importlib.util.spec_from_file_location("sascorer", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load sascorer from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules["sascorer"] = module
    spec.loader.exec_module(module)
    _SASCORER = module
    return module


def sa_score(mol: Chem.Mol) -> float:
    """Return the real Ertl synthetic accessibility score for *mol*."""
    return float(load_sascorer().calculateScore(mol))


def compute_properties(mol: Chem.Mol, config: PropertiesConfig | None = None) -> dict[str, Any]:
    """Compute descriptors and drug-likeness flags for a sanitized molecule."""
    config = config or PropertiesConfig()
    mw = float(Descriptors.MolWt(mol))
    clogp = float(Crippen.MolLogP(mol))
    tpsa = float(rdMolDescriptors.CalcTPSA(mol))
    hbd = int(rdMolDescriptors.CalcNumLipinskiHBD(mol))
    hba = int(rdMolDescriptors.CalcNumLipinskiHBA(mol))
    rotatable = int(rdMolDescriptors.CalcNumRotatableBonds(mol))
    fraction_csp3 = float(rdMolDescriptors.CalcFractionCSP3(mol))
    ring_count = int(rdMolDescriptors.CalcNumRings(mol))
    qed = float(QED.qed(mol))
    sa = sa_score(mol)

    violations = 0
    if mw > config.mw_limit:
        violations += 1
    if clogp > config.clogp_limit:
        violations += 1
    if hbd > config.hbd_limit:
        violations += 1
    if hba > config.hba_limit:
        violations += 1
    lipinski_pass = violations <= config.lipinski_pass_max_violations
    lipinski_strict_pass = violations == 0
    veber_pass = (
        rotatable <= config.veber_rotatable_bonds_limit and tpsa <= config.veber_tpsa_limit
    )
    return {
        "mw": mw,
        "clogp": clogp,
        "tpsa": tpsa,
        "hbd": hbd,
        "hba": hba,
        "rotatable_bonds": rotatable,
        "fraction_csp3": fraction_csp3,
        "ring_count": ring_count,
        "qed": qed,
        "sa_score": sa,
        "lipinski_violations": violations,
        "lipinski_pass": lipinski_pass,
        "lipinski_strict_pass": lipinski_strict_pass,
        "veber_pass": veber_pass,
    }


def property_record(
    compound_id: str,
    smiles: str,
    config: PropertiesConfig | None = None,
    *,
    label: str | None = None,
) -> PropertyRecord:
    config = config or PropertiesConfig()
    raw = "" if smiles is None else str(smiles).strip()
    if not raw:
        return PropertyRecord(
            compound_id=str(compound_id),
            smiles="",
            label=label,
            error="empty SMILES",
        )
    mol = Chem.MolFromSmiles(raw)
    if mol is None:
        return PropertyRecord(
            compound_id=str(compound_id),
            smiles=raw,
            label=label,
            error="invalid SMILES",
        )
    try:
        values = compute_properties(mol, config)
    except Exception as exc:  # noqa: BLE001
        return PropertyRecord(
            compound_id=str(compound_id),
            smiles=raw,
            label=label,
            error=f"{type(exc).__name__}: {exc}",
        )
    return PropertyRecord(compound_id=str(compound_id), smiles=raw, label=label, **values)


def _median(values: list[float]) -> float | None:
    return float(median(values)) if values else None


def _label_value(label: str | None) -> str | None:
    if label is None:
        return None
    text = str(label).strip()
    return text or None


def summarize_properties(records: list[PropertyRecord]) -> PropertiesSummary:
    valid = [rec for rec in records if rec.error is None and rec.mw is not None]
    def col(name: str) -> list[float]:
        return [float(getattr(rec, name)) for rec in valid if getattr(rec, name) is not None]

    n_valid = len(valid)
    lip_pass = sum(1 for rec in valid if rec.lipinski_pass)
    lip_strict = sum(1 for rec in valid if rec.lipinski_strict_pass)
    veber = sum(1 for rec in valid if rec.veber_pass)

    by_label: dict[str, PropertyLabelStats] = {}
    grouped: dict[str, list[PropertyRecord]] = defaultdict(list)
    for rec in records:
        key = _label_value(rec.label)
        if key:
            grouped[key].append(rec)
    for label, subset in sorted(grouped.items()):
        sub_valid = [rec for rec in subset if rec.error is None and rec.mw is not None]
        n_sub = len(sub_valid)
        sub_lip = sum(1 for rec in sub_valid if rec.lipinski_pass)
        sub_strict = sum(1 for rec in sub_valid if rec.lipinski_strict_pass)
        sub_veber = sum(1 for rec in sub_valid if rec.veber_pass)
        by_label[label] = PropertyLabelStats(
            n_total=len(subset),
            n_valid=n_sub,
            median_mw=_median([float(rec.mw) for rec in sub_valid if rec.mw is not None]),
            median_clogp=_median([float(rec.clogp) for rec in sub_valid if rec.clogp is not None]),
            median_tpsa=_median([float(rec.tpsa) for rec in sub_valid if rec.tpsa is not None]),
            median_qed=_median([float(rec.qed) for rec in sub_valid if rec.qed is not None]),
            median_sa_score=_median(
                [float(rec.sa_score) for rec in sub_valid if rec.sa_score is not None]
            ),
            lipinski_pass_count=sub_lip,
            lipinski_pass_rate=(sub_lip / n_sub) if n_sub else 0.0,
            lipinski_strict_pass_count=sub_strict,
            lipinski_strict_pass_rate=(sub_strict / n_sub) if n_sub else 0.0,
            veber_pass_count=sub_veber,
            veber_pass_rate=(sub_veber / n_sub) if n_sub else 0.0,
        )

    return PropertiesSummary(
        n_total=len(records),
        n_valid=n_valid,
        n_invalid=len(records) - n_valid,
        median_mw=_median(col("mw")),
        median_clogp=_median(col("clogp")),
        median_tpsa=_median(col("tpsa")),
        median_hbd=_median(col("hbd")),
        median_hba=_median(col("hba")),
        median_rotatable_bonds=_median(col("rotatable_bonds")),
        median_fraction_csp3=_median(col("fraction_csp3")),
        median_ring_count=_median(col("ring_count")),
        median_qed=_median(col("qed")),
        median_sa_score=_median(col("sa_score")),
        lipinski_pass_count=lip_pass,
        lipinski_pass_rate=(lip_pass / n_valid) if n_valid else 0.0,
        lipinski_strict_pass_count=lip_strict,
        lipinski_strict_pass_rate=(lip_strict / n_valid) if n_valid else 0.0,
        veber_pass_count=veber,
        veber_pass_rate=(veber / n_valid) if n_valid else 0.0,
        by_label=by_label,
    )


def run_properties(
    rows: Iterable[dict[str, Any]],
    config: PropertiesConfig | None = None,
    *,
    smiles_field: str = "smiles",
    id_field: str = "compound_id",
    label_field: str = "label",
) -> tuple[list[PropertyRecord], PropertiesSummary]:
    config = config or PropertiesConfig()
    records: list[PropertyRecord] = []
    for index, row in enumerate(rows):
        compound_id = row.get(id_field) or row.get("compound_id") or f"cmpd_{index + 1:06d}"
        smiles = row.get(smiles_field, row.get("smiles", ""))
        raw_label = row.get(label_field, row.get("label"))
        label = _label_value(None if raw_label == "" else raw_label)
        records.append(property_record(str(compound_id), str(smiles or ""), config, label=label))
    return records, summarize_properties(records)


def run_properties_file(
    input_path: str | Path,
    output_dir: str | Path,
    config: PropertiesConfig | None = None,
    *,
    smiles_column: str | None = None,
    id_column: str | None = None,
    label_column: str | None = None,
):
    """Compute properties for a compound table and write CSV plus audit log."""
    config = config or PropertiesConfig()
    input_path = Path(input_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    frame = load_compound_frame(
        input_path,
        smiles_column=smiles_column,
        id_column=id_column,
        label_column=label_column,
    )
    records, summary = run_properties(frame.to_dict(orient="records"), config)
    csv_path = output_dir / "properties.csv"
    write_csv(csv_path, [rec.to_csv_row() for rec in records])
    manifest = build_audit(
        command="properties",
        config=config,
        input_path=input_path,
        output_dir=output_dir,
        summary=summary,
        outputs=[str(csv_path.resolve()), str((output_dir / "audit_log.json").resolve())],
    )
    write_audit(output_dir / "audit_log.json", manifest)
    return records, summary, manifest
