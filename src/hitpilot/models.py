"""Pydantic v2 models for HitPilot hygiene, properties, novelty, and diversity."""

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



# --- Phase 2: properties, novelty, and diversity ---------------------------------

class PropertiesConfig(BaseModel):
    """Thresholds for Lipinski and Veber drug-likeness flags."""

    mw_limit: float = 500.0
    clogp_limit: float = 5.0
    hbd_limit: int = 5
    hba_limit: int = 10
    lipinski_pass_max_violations: int = 1
    veber_rotatable_bonds_limit: int = 10
    veber_tpsa_limit: float = 140.0


class PropertyRecord(BaseModel):
    """Per-compound RDKit physicochemical descriptors and drug-likeness flags."""

    compound_id: str
    smiles: str
    mw: float | None = None
    clogp: float | None = None
    tpsa: float | None = None
    hbd: int | None = None
    hba: int | None = None
    rotatable_bonds: int | None = None
    fraction_csp3: float | None = None
    ring_count: int | None = None
    qed: float | None = None
    sa_score: float | None = None
    lipinski_violations: int | None = None
    lipinski_pass: bool | None = None
    lipinski_strict_pass: bool | None = None
    veber_pass: bool | None = None
    label: str | None = None
    error: str | None = None

    def to_csv_row(self) -> dict[str, object]:
        return {
            "compound_id": self.compound_id,
            "smiles": self.smiles,
            "mw": "" if self.mw is None else self.mw,
            "clogp": "" if self.clogp is None else self.clogp,
            "tpsa": "" if self.tpsa is None else self.tpsa,
            "hbd": "" if self.hbd is None else self.hbd,
            "hba": "" if self.hba is None else self.hba,
            "rotatable_bonds": "" if self.rotatable_bonds is None else self.rotatable_bonds,
            "fraction_csp3": "" if self.fraction_csp3 is None else self.fraction_csp3,
            "ring_count": "" if self.ring_count is None else self.ring_count,
            "qed": "" if self.qed is None else self.qed,
            "sa_score": "" if self.sa_score is None else self.sa_score,
            "lipinski_violations": (
                "" if self.lipinski_violations is None else self.lipinski_violations
            ),
            "lipinski_pass": "" if self.lipinski_pass is None else self.lipinski_pass,
            "lipinski_strict_pass": (
                "" if self.lipinski_strict_pass is None else self.lipinski_strict_pass
            ),
            "veber_pass": "" if self.veber_pass is None else self.veber_pass,
            "label": self.label or "",
            "error": self.error or "",
        }


class PropertyLabelStats(BaseModel):
    """Property aggregates restricted to a single activity label."""

    n_total: int = 0
    n_valid: int = 0
    median_mw: float | None = None
    median_clogp: float | None = None
    median_tpsa: float | None = None
    median_qed: float | None = None
    median_sa_score: float | None = None
    lipinski_pass_count: int = 0
    lipinski_pass_rate: float = 0.0
    lipinski_strict_pass_count: int = 0
    lipinski_strict_pass_rate: float = 0.0
    veber_pass_count: int = 0
    veber_pass_rate: float = 0.0


class PropertiesSummary(BaseModel):
    """Aggregate physicochemical statistics for a properties run."""

    n_total: int
    n_valid: int
    n_invalid: int
    median_mw: float | None = None
    median_clogp: float | None = None
    median_tpsa: float | None = None
    median_hbd: float | None = None
    median_hba: float | None = None
    median_rotatable_bonds: float | None = None
    median_fraction_csp3: float | None = None
    median_ring_count: float | None = None
    median_qed: float | None = None
    median_sa_score: float | None = None
    lipinski_pass_count: int = 0
    lipinski_pass_rate: float = 0.0
    lipinski_strict_pass_count: int = 0
    lipinski_strict_pass_rate: float = 0.0
    veber_pass_count: int = 0
    veber_pass_rate: float = 0.0
    by_label: dict[str, PropertyLabelStats] = Field(default_factory=dict)


class ChemblActive(BaseModel):
    """One unique ChEMBL reference bioactive after SMILES deduplication."""

    molecule_chembl_id: str
    smiles: str
    canonical_smiles: str
    pchembl_value: float
    target_chembl_id: str


class NoveltyConfig(BaseModel):
    """Morgan-fingerprint novelty check against ChEMBL reference actives."""

    target_chembl_id: str = "CHEMBL206"
    pchembl_min: float = 6.0
    max_records: int = 2000
    cache_dir: str = "data/chembl"
    refresh_cache: bool = False
    radius: int = 2
    fp_size: int = 2048
    known_match_min: float = 0.85
    close_analog_min: float = 0.65
    moderate_novelty_min: float = 0.45
    timeout: int = 15


class NoveltyRecord(BaseModel):
    """Nearest-neighbor novelty result for one compound."""

    compound_id: str
    smiles: str
    max_tanimoto: float | None = None
    nearest_ref_id: str | None = None
    nearest_ref_pchembl: float | None = None
    nearest_ref_smiles: str | None = None
    novelty_tier: str | None = None
    label: str | None = None
    error: str | None = None

    def to_csv_row(self) -> dict[str, object]:
        return {
            "compound_id": self.compound_id,
            "smiles": self.smiles,
            "max_tanimoto": "" if self.max_tanimoto is None else self.max_tanimoto,
            "nearest_ref_id": self.nearest_ref_id or "",
            "nearest_ref_pchembl": (
                "" if self.nearest_ref_pchembl is None else self.nearest_ref_pchembl
            ),
            "nearest_ref_smiles": self.nearest_ref_smiles or "",
            "novelty_tier": self.novelty_tier or "",
            "label": self.label or "",
            "error": self.error or "",
        }


class NoveltySummary(BaseModel):
    """Aggregate novelty-tier counts."""

    n_total: int
    n_valid: int
    n_invalid: int
    n_reference_actives: int
    target_chembl_id: str
    pchembl_min: float
    median_max_tanimoto: float | None = None
    counts_by_tier: dict[str, int] = Field(default_factory=dict)
    by_label: dict[str, dict[str, int]] = Field(default_factory=dict)


class PubChemLookupResult(BaseModel):
    """Result of a targeted PubChem PUG identity lookup."""

    query_smiles: str
    cids: list[int] = Field(default_factory=list)
    error: str | None = None


class PubChemSimilarResult(BaseModel):
    """Result of a targeted PubChem PUG 2D similarity search."""

    query_smiles: str
    threshold: int
    max_records: int
    cids: list[int] = Field(default_factory=list)
    error: str | None = None


class DiversityConfig(BaseModel):
    """Scaffold, Butina clustering, and sphere-exclusion settings."""

    radius: int = 2
    fp_size: int = 2048
    tanimoto_cutoff: float = 0.60
    exclusion_similarity: float = 0.65
    score_column: str | None = None
    score_ascending: bool = False
    rank_by_score: bool = False


class DiversityRecord(BaseModel):
    """Scaffold, cluster, and sphere-exclusion assignment for one compound."""

    compound_id: str
    smiles: str
    murcko_scaffold_smiles: str = ""
    generic_scaffold_smiles: str = ""
    acyclic: bool = False
    scaffold_frequency: int = 0
    cluster_id: int | None = None
    cluster_size: int = 0
    is_cluster_centroid: bool = False
    sphere_selected: bool = False
    sphere_rank: int | None = None
    excluded_by: str | None = None
    label: str | None = None
    error: str | None = None

    def to_csv_row(self) -> dict[str, object]:
        return {
            "compound_id": self.compound_id,
            "smiles": self.smiles,
            "murcko_scaffold_smiles": self.murcko_scaffold_smiles,
            "generic_scaffold_smiles": self.generic_scaffold_smiles,
            "acyclic": self.acyclic,
            "scaffold_frequency": self.scaffold_frequency,
            "cluster_id": "" if self.cluster_id is None else self.cluster_id,
            "cluster_size": self.cluster_size,
            "is_cluster_centroid": self.is_cluster_centroid,
            "sphere_selected": self.sphere_selected,
            "sphere_rank": "" if self.sphere_rank is None else self.sphere_rank,
            "excluded_by": self.excluded_by or "",
            "label": self.label or "",
            "error": self.error or "",
        }


class DiversitySummary(BaseModel):
    """Scaffold, cluster, and diverse-subset counts."""

    n_total: int
    n_valid: int
    n_invalid: int
    n_unique_murcko_scaffolds: int = 0
    n_unique_generic_scaffolds: int = 0
    n_acyclic: int = 0
    n_clusters: int = 0
    top_cluster_sizes: list[int] = Field(default_factory=list)
    n_sphere_selected: int = 0
    n_actives_selected: int | None = None
    n_actives_total: int | None = None
    tanimoto_cutoff: float
    exclusion_similarity: float


class PhaseAuditManifest(BaseModel):
    """Machine-readable audit trail for a Phase 2 command."""

    timestamp_utc: str
    hitpilot_version: str
    rdkit_version: str
    command: str
    config: dict = Field(default_factory=dict)
    input_path: str
    input_sha256: str
    output_dir: str
    summary: dict = Field(default_factory=dict)
    outputs: list[str] = Field(default_factory=list)
