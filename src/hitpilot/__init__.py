"""HitPilot: auditable hit-triage and med-chem copilot."""

from hitpilot.benchmark import run_litpcba_benchmark
from hitpilot.diversity import run_diversity, run_diversity_file
from hitpilot.dossier import (
    DossierConfig,
    DossierManifest,
    DossierResult,
    generate_medchem_dossier,
    generate_medchem_dossier_from_dirs,
)
from hitpilot.hygiene import run_hygiene, run_hygiene_file, standardize_smiles
from hitpilot.mcp_server import create_chem_mcp_server
from hitpilot.models import (
    AlertHit,
    AuditManifest,
    DiversityConfig,
    DiversityRecord,
    DiversitySummary,
    HygieneConfig,
    HygieneRecord,
    HygieneSummary,
    LabelBreakdown,
    NoveltyConfig,
    NoveltyRecord,
    NoveltySummary,
    PropertiesConfig,
    PropertiesSummary,
    PropertyRecord,
)
from hitpilot.novelty import (
    fetch_chembl_target_actives,
    lookup_pubchem_compound,
    run_novelty,
    run_novelty_file,
    search_pubchem_similar,
)
from hitpilot.phase2 import run_phase2_file
from hitpilot.properties import run_properties, run_properties_file, sa_score
from hitpilot.scoring import (
    EnrichmentMetrics,
    ScoredShortlistRecord,
    ScorerHookConfig,
    ScoringSummary,
    compute_enrichment_metrics,
    generate_litpcba_baseline_scores,
    merge_external_scores,
    run_scorer_hook,
    run_scorer_hook_file,
    select_diverse_shortlist,
)
from hitpilot.version import __version__

__all__ = [
    "__version__",
    "AlertHit",
    "AuditManifest",
    "DiversityConfig",
    "DiversityRecord",
    "DiversitySummary",
    "EnrichmentMetrics",
    "HygieneConfig",
    "HygieneRecord",
    "HygieneSummary",
    "LabelBreakdown",
    "NoveltyConfig",
    "NoveltyRecord",
    "NoveltySummary",
    "PropertiesConfig",
    "PropertyRecord",
    "PropertiesSummary",
    "ScoredShortlistRecord",
    "ScorerHookConfig",
    "ScoringSummary",
    "DossierConfig",
    "DossierManifest",
    "DossierResult",
    "compute_enrichment_metrics",
    "create_chem_mcp_server",
    "fetch_chembl_target_actives",
    "generate_litpcba_baseline_scores",
    "generate_medchem_dossier",
    "generate_medchem_dossier_from_dirs",
    "lookup_pubchem_compound",
    "merge_external_scores",
    "run_diversity",
    "run_diversity_file",
    "run_hygiene",
    "run_hygiene_file",
    "run_litpcba_benchmark",
    "run_novelty",
    "run_novelty_file",
    "run_phase2_file",
    "run_properties",
    "run_properties_file",
    "run_scorer_hook",
    "run_scorer_hook_file",
    "sa_score",
    "search_pubchem_similar",
    "select_diverse_shortlist",
    "standardize_smiles",
]
