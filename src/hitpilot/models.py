"""Pydantic v2 models for HitPilot Phase 1 hygiene."""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, Field

DedupKey = Literal["canonical_smiles", "inchikey"]


class HygieneConfig(BaseModel):
    """Controls standardization, deduplication, and FilterCatalog behavior.

    PAINS, Brenk, and NIH alerts are enabled and dropped by default. Set the
    matching ``drop_*`` flag to ``False`` to record the alert without removing
    the compound (flag-only mode).
    """

    canonicalize_tautomers: bool = False
    dedup_key: DedupKey = "canonical_smiles"
    score_column: str | None = None
    score_ascending: bool = False
    enable_pains: bool = True
    enable_brenk: bool = True
    enable_nih: bool = True
    drop_pains: bool = True
    drop_brenk: bool = True
    drop_nih: bool = True


class AlertHit(BaseModel):
    """One RDKit FilterCatalog match on a compound."""

    catalog: str = Field(description="Normalized catalog name, e.g. PAINS_A or BRENK.")
    rule: str = Field(description="Rule description from the FilterCatalog entry.")


class HygieneRecord(BaseModel):
    """Per-compound hygiene result."""

    compound_id: str
    input_smiles: str
    standardized_smiles: str | None = None
    inchikey: str | None = None
    passed: bool = False
    drop_reasons: list[str] = Field(default_factory=list)
    alerts: list[AlertHit] = Field(default_factory=list)
    duplicate_of: str | None = None
    label: str | None = None
    score: float | None = None
    error: str | None = None
    input_index: int = 0

    def alerts_text(self) -> str:
        return "|".join(f"{hit.catalog}:{hit.rule}" for hit in self.alerts)

    def drop_reasons_text(self) -> str:
        return "|".join(self.drop_reasons)

    def to_csv_row(self) -> dict[str, object]:
        return {
            "compound_id": self.compound_id,
            "input_smiles": self.input_smiles,
            "standardized_smiles": self.standardized_smiles or "",
            "inchikey": self.inchikey or "",
            "passed": self.passed,
            "drop_reasons": self.drop_reasons_text(),
            "alerts": self.alerts_text(),
            "duplicate_of": self.duplicate_of or "",
            "label": self.label or "",
            "score": "" if self.score is None else self.score,
            "error": self.error or "",
        }


class LabelBreakdown(BaseModel):
    """Hygiene counts restricted to a single activity label."""

    n_total: int = 0
    n_passed: int = 0
    n_dropped: int = 0
    n_with_alerts: int = 0
    pass_rate: float = 0.0
    counts_by_drop_reason: dict[str, int] = Field(default_factory=dict)


class HygieneSummary(BaseModel):
    """Aggregate counts for a hygiene run."""

    n_total: int
    n_passed: int
    n_dropped: int
    pass_rate: float
    counts_by_drop_reason: dict[str, int] = Field(default_factory=dict)
    counts_by_catalog: dict[str, int] = Field(default_factory=dict)
    by_label: dict[str, LabelBreakdown] = Field(default_factory=dict)


class AuditManifest(BaseModel):
    """Machine-readable audit trail for one hygiene run."""

    timestamp_utc: str
    hitpilot_version: str
    rdkit_version: str
    config: HygieneConfig
    input_path: str
    input_sha256: str
    output_dir: str
    summary: HygieneSummary
