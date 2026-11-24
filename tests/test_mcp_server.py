"""Tests for ChemMCP tool registration and in-process tool invocation."""

from __future__ import annotations

import asyncio
import json
from pathlib import Path

from mcp import Client

from hitpilot.mcp_server import create_chem_mcp_server
from hitpilot.models import ChemblActive
from hitpilot.novelty import canonicalize_smiles

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
IBUPROFEN = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
BENZENE = "c1ccccc1"
ETHANOL = "CCO"
EXPECTED_TOOLS = {
    "standardize_and_screen_smiles",
    "compute_physchem_properties",
    "check_chembl_novelty",
    "lookup_pubchem_compound",
    "analyze_scaffolds_and_diversity",
    "evaluate_rescorer_csv",
    "run_hygiene_csv",
    "run_phase2_csv",
    "generate_medchem_dossier",
}


def _payload(result) -> dict:
    if getattr(result, "is_error", False):
        raise AssertionError(f"tool error: {result}")
    if getattr(result, "structured_content", None):
        return result.structured_content
    text = result.content[0].text
    return json.loads(text)


def test_chemmcp_registers_nine_tools_with_schemas():
    async def _run():
        server = create_chem_mcp_server()
        assert server.name == "ChemMCP"
        async with Client(server) as client:
            listing = await client.list_tools()
            tools = listing.tools
            names = {tool.name for tool in tools}
            assert names == EXPECTED_TOOLS
            for tool in tools:
                schema = tool.input_schema
                assert isinstance(schema, dict)
                assert schema.get("type") == "object"
                assert "properties" in schema
                assert isinstance(schema["properties"], dict)
                assert schema["properties"], tool.name

    asyncio.run(_run())


def test_chemmcp_standardize_and_properties_tools():
    async def _run():
        server = create_chem_mcp_server()
        async with Client(server) as client:
            screened = _payload(
                await client.call_tool(
                    "standardize_and_screen_smiles",
                    {"smiles_list": [IBUPROFEN, "CCNCC.Cl", "not_a_smiles"]},
                )
            )
            assert "rdkit_version" in screened
            recs = {row["compound_id"]: row for row in screened["records"]}
            assert recs["smi_0001"]["passed"] is True
            assert recs["smi_0002"]["standardized_smiles"]
            assert "Cl" not in (recs["smi_0002"]["standardized_smiles"] or "")
            assert recs["smi_0003"]["drop_reasons"] == ["invalid_smiles"]
            assert screened["summary"]["n_total"] == 3

            props = _payload(
                await client.call_tool(
                    "compute_physchem_properties",
                    {"smiles_list": [ASPIRIN]},
                )
            )
            rec = props["records"][0]
            assert abs(rec["mw"] - 180.159) < 0.01
            assert rec["lipinski_strict_pass"] is True
            assert rec["veber_pass"] is True
            assert rec["hbd"] == 1
            assert rec["sa_score"] > 1

    asyncio.run(_run())


def test_chemmcp_diversity_and_evaluate_rescorer(tmp_path: Path):
    async def _run():
        server = create_chem_mcp_server()
        async with Client(server) as client:
            diversity = _payload(
                await client.call_tool(
                    "analyze_scaffolds_and_diversity",
                    {
                        "smiles_list": [BENZENE, "C1=CC=CC=C1", IBUPROFEN, ETHANOL],
                        "butina_cutoff": 0.60,
                        "exclusion_similarity": 0.65,
                    },
                )
            )
            assert diversity["summary"]["n_total"] == 4
            assert diversity["summary"]["n_acyclic"] == 1
            selected = [row for row in diversity["records"] if row["sphere_selected"]]
            ids = {row["compound_id"] for row in selected}
            # The first benzene should be selected; the kekule form should not.
            assert "smi_0001" in ids
            assert "smi_0002" not in ids

            compounds = tmp_path / "compounds.csv"
            scores = tmp_path / "scores.csv"
            compounds.write_text(
                "compound_id,smiles,label\n"
                f"ibu,{IBUPROFEN},active\n"
                f"asp,{ASPIRIN},inactive\n",
                encoding="utf-8",
            )
            scores.write_text("compound_id,score\nibu,0.9\nasp,0.1\n", encoding="utf-8")
            out = tmp_path / "hook"
            evaluated = _payload(
                await client.call_tool(
                    "evaluate_rescorer_csv",
                    {
                        "compounds_csv": str(compounds),
                        "scores_csv": str(scores),
                        "output_dir": str(out),
                        "higher_is_better": True,
                        "top_n": 2,
                    },
                )
            )
            assert evaluated["summary"]["enrichment"]["roc_auc"] == 1.0
            assert (out / "ranked_compounds.csv").is_file()
            assert evaluated["outputs"]

            hygiene = _payload(
                await client.call_tool(
                    "run_hygiene_csv",
                    {"input_csv": str(compounds), "output_dir": str(tmp_path / "hyg")},
                )
            )
            assert hygiene["summary"]["n_total"] == 2
            assert (tmp_path / "hyg" / "passed_compounds.csv").is_file()

    asyncio.run(_run())


def test_chemmcp_novelty_and_pubchem_are_mocked(monkeypatch):
    from hitpilot import mcp_server as mcp_mod
    from hitpilot.models import PubChemLookupResult, PubChemSimilarResult

    refs = [
        ChemblActive(
            molecule_chembl_id="CHEMBL521",
            smiles=IBUPROFEN,
            canonical_smiles=canonicalize_smiles(IBUPROFEN),
            pchembl_value=7.2,
            target_chembl_id="CHEMBL206",
        )
    ]
    monkeypatch.setattr(mcp_mod, "fetch_chembl_target_actives", lambda **kwargs: refs)
    monkeypatch.setattr(
        mcp_mod,
        "_lookup_pubchem_compound",
        lambda smiles, **kwargs: PubChemLookupResult(query_smiles=smiles, cids=[2244]),
    )
    monkeypatch.setattr(
        mcp_mod,
        "_search_pubchem_similar",
        lambda smiles, threshold=85, max_records=5, **kwargs: PubChemSimilarResult(
            query_smiles=smiles,
            threshold=threshold,
            max_records=max_records,
            cids=[2244, 3672],
        ),
    )

    async def _run():
        server = create_chem_mcp_server()
        async with Client(server) as client:
            novelty = _payload(
                await client.call_tool(
                    "check_chembl_novelty",
                    {"smiles_list": [IBUPROFEN], "target_chembl_id": "CHEMBL206"},
                )
            )
            rec = novelty["records"][0]
            assert rec["nearest_ref_id"] == "CHEMBL521"
            assert rec["max_tanimoto"] == 1.0
            assert rec["novelty_tier"] == "known_or_close_match"

            pubchem = _payload(
                await client.call_tool(
                    "lookup_pubchem_compound",
                    {"smiles": ASPIRIN, "similarity_threshold": 85, "max_similar": 5},
                )
            )
            assert pubchem["cids"] == [2244]
            assert pubchem["similar_cids"] == [2244, 3672]

    asyncio.run(_run())
