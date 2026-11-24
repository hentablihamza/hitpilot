"""ChemMCP: stdio MCP server exposing verified HitPilot hit-triage tools."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from mcp.server.mcpserver import MCPServer

from hitpilot.diversity import run_diversity
from hitpilot.hygiene import run_hygiene, run_hygiene_file
from hitpilot.io_utils import rdkit_version
from hitpilot.models import DiversityConfig, HygieneConfig, NoveltyConfig
from hitpilot.novelty import (
    fetch_chembl_target_actives,
    run_novelty,
)
from hitpilot.novelty import (
    lookup_pubchem_compound as _lookup_pubchem_compound,
)
from hitpilot.novelty import (
    search_pubchem_similar as _search_pubchem_similar,
)
from hitpilot.phase2 import run_phase2_file
from hitpilot.properties import run_properties
from hitpilot.scoring import ScorerHookConfig, run_scorer_hook_file
from hitpilot.version import __version__

INSTRUCTIONS = (
    "ChemMCP provides verifiable RDKit, ChEMBL, and PubChem hit-triage tools "
    "from HitPilot. All chemistry is computed with real RDKit code or live "
    "public APIs; tools never invent SMILES, scores, or metrics."
)


def _rows_from_smiles(smiles_list: list[str]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for index, smiles in enumerate(smiles_list, start=1):
        rows.append({"compound_id": f"smi_{index:04d}", "smiles": str(smiles)})
    return rows


def create_chem_mcp_server() -> MCPServer:
    """Build a ChemMCP server with all eight typed hit-triage tools registered."""
    server = MCPServer(
        name="ChemMCP",
        instructions=INSTRUCTIONS,
        version=__version__,
    )

    @server.tool()
    def standardize_and_screen_smiles(
        smiles_list: list[str],
        drop_pains: bool = True,
        drop_brenk: bool = True,
        drop_nih: bool = True,
    ) -> dict[str, Any]:
        """Sanitize, strip salts, neutralize, canonicalize, deduplicate, and screen alerts."""
        config = HygieneConfig(
            drop_pains=drop_pains,
            drop_brenk=drop_brenk,
            drop_nih=drop_nih,
        )
        records, summary = run_hygiene(_rows_from_smiles(smiles_list), config)
        return {
            "rdkit_version": rdkit_version(),
            "hitpilot_version": __version__,
            "records": [rec.model_dump(mode="json") for rec in records],
            "summary": summary.model_dump(mode="json"),
        }

    @server.tool()
    def compute_physchem_properties(smiles_list: list[str]) -> dict[str, Any]:
        """Compute RDKit descriptors, Ertl SA score, and Lipinski/Veber flags."""
        records, summary = run_properties(_rows_from_smiles(smiles_list))
        return {
            "rdkit_version": rdkit_version(),
            "hitpilot_version": __version__,
            "records": [rec.model_dump(mode="json") for rec in records],
            "summary": summary.model_dump(mode="json"),
        }

    @server.tool()
    def check_chembl_novelty(
        smiles_list: list[str],
        target_chembl_id: str = "CHEMBL206",
        pchembl_min: float = 6.0,
        max_records: int = 2000,
    ) -> dict[str, Any]:
        """Morgan ECFP4 nearest-neighbor novelty versus ChEMBL target actives."""
        config = NoveltyConfig(
            target_chembl_id=target_chembl_id,
            pchembl_min=pchembl_min,
            max_records=max_records,
        )
        references = fetch_chembl_target_actives(
            target_chembl_id=target_chembl_id,
            pchembl_min=pchembl_min,
            max_records=max_records,
            cache_dir=config.cache_dir,
        )
        records, summary = run_novelty(_rows_from_smiles(smiles_list), references, config)
        return {
            "rdkit_version": rdkit_version(),
            "hitpilot_version": __version__,
            "records": [rec.model_dump(mode="json") for rec in records],
            "summary": summary.model_dump(mode="json"),
        }

    @server.tool()
    def lookup_pubchem_compound(
        smiles: str,
        similarity_threshold: int = 85,
        max_similar: int = 5,
    ) -> dict[str, Any]:
        """Look up a PubChem CID and 2D similar CIDs via PUG REST."""
        identity = _lookup_pubchem_compound(smiles)
        similar = _search_pubchem_similar(
            smiles,
            threshold=similarity_threshold,
            max_records=max_similar,
        )
        return {
            "query_smiles": smiles,
            "cids": identity.cids,
            "identity_error": identity.error,
            "similar_cids": similar.cids,
            "similar_error": similar.error,
            "similarity_threshold": similarity_threshold,
            "max_similar": max_similar,
        }

    @server.tool()
    def analyze_scaffolds_and_diversity(
        smiles_list: list[str],
        butina_cutoff: float = 0.60,
        exclusion_similarity: float = 0.65,
    ) -> dict[str, Any]:
        """Return Bemis-Murcko scaffolds, Butina clusters, and a diverse subset."""
        config = DiversityConfig(
            tanimoto_cutoff=butina_cutoff,
            exclusion_similarity=exclusion_similarity,
        )
        records, summary = run_diversity(_rows_from_smiles(smiles_list), config)
        return {
            "rdkit_version": rdkit_version(),
            "hitpilot_version": __version__,
            "records": [rec.model_dump(mode="json") for rec in records],
            "summary": summary.model_dump(mode="json"),
        }

    @server.tool()
    def evaluate_rescorer_csv(
        compounds_csv: str,
        scores_csv: str,
        output_dir: str,
        score_column: str = "score",
        id_column: str = "compound_id",
        higher_is_better: bool = True,
        top_n: int = 25,
        exclusion_similarity: float = 0.65,
    ) -> dict[str, Any]:
        """Run the pluggable scorer hook and return enrichment plus shortlist preview."""
        config = ScorerHookConfig(
            id_column=id_column,
            score_column=score_column,
            higher_is_better=higher_is_better,
            top_n=top_n,
            diversity_cutoff=exclusion_similarity,
        )
        _ranked, shortlist, summary, manifest = run_scorer_hook_file(
            compounds_csv,
            scores_csv,
            output_dir,
            config,
        )
        preview = shortlist.head(min(top_n, 10)).to_dict(orient="records")
        return {
            "rdkit_version": rdkit_version(),
            "hitpilot_version": __version__,
            "summary": summary.model_dump(mode="json"),
            "shortlist_preview": preview,
            "outputs": manifest.outputs,
            "output_dir": str(Path(output_dir).resolve()),
        }

    @server.tool()
    def run_hygiene_csv(
        input_csv: str,
        output_dir: str,
        drop_pains: bool = True,
        drop_brenk: bool = True,
        drop_nih: bool = True,
    ) -> dict[str, Any]:
        """Run the Phase 1 hygiene batch pipeline on a CSV/TSV/SMI file."""
        config = HygieneConfig(
            drop_pains=drop_pains,
            drop_brenk=drop_brenk,
            drop_nih=drop_nih,
        )
        manifest = run_hygiene_file(input_csv, output_dir, config)
        return {
            "rdkit_version": rdkit_version(),
            "hitpilot_version": __version__,
            "summary": manifest.summary.model_dump(mode="json"),
            "output_dir": str(Path(output_dir).resolve()),
            "outputs": [
                str(Path(output_dir) / "passed_compounds.csv"),
                str(Path(output_dir) / "all_compounds.csv"),
                str(Path(output_dir) / "audit_log.json"),
            ],
            "input_sha256": manifest.input_sha256,
        }

    @server.tool()
    def run_phase2_csv(
        input_csv: str,
        output_dir: str,
        target_chembl_id: str = "CHEMBL206",
        pchembl_min: float = 6.0,
        max_records: int = 2000,
    ) -> dict[str, Any]:
        """Run the Phase 2 properties, novelty, and diversity pipeline on a file."""
        result = run_phase2_file(
            input_csv,
            output_dir,
            novelty_config=NoveltyConfig(
                target_chembl_id=target_chembl_id,
                pchembl_min=pchembl_min,
                max_records=max_records,
            ),
        )
        manifest = result["manifest"]
        return {
            "rdkit_version": rdkit_version(),
            "hitpilot_version": __version__,
            "summary": manifest.summary,
            "output_dir": str(Path(output_dir).resolve()),
            "outputs": manifest.outputs,
            "input_sha256": manifest.input_sha256,
        }

    return server


def main() -> None:
    """Run ChemMCP over stdio."""
    create_chem_mcp_server().run(transport="stdio")


if __name__ == "__main__":
    main()
