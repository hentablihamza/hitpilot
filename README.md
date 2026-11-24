# HitPilot

HitPilot is an open-source, auditable hit-triage and medicinal-chemistry copilot.
Phase 1 standardizes public screening libraries with RDKit, strips salts, drops
exact duplicates, and flags or removes PAINS, Brenk, and NIH interference
alerts — with a machine-readable audit log for every compound.

The public benchmark target is a LIT-PCBA protein (real PubChem BioAssay
actives and inactives, not property-matched decoys).

## Phase 1 quickstart

All commands use the project virtual environment:

```bash
# Editable install (dependencies already present in .venv)
uv pip install --python .venv/bin/python -e . --no-deps

# Fetch one public LIT-PCBA target (stores files under data/litpcba/)
.venv/bin/hitpilot fetch-litpcba ESR1_ant

# Run hygiene on the unified actives + inactives table
.venv/bin/hitpilot hygiene data/litpcba/ESR1_ant/esr1_ant_litpcba.csv \
  --output-dir outputs/phase1_ESR1_ant

# Tests and lint
.venv/bin/pytest -v
.venv/bin/ruff check src tests
```

Outputs written to the chosen directory:

- `passed_compounds.csv` — compounds that survived standardization, dedup, and drop filters
- `all_compounds.csv` — every input row with standardized SMILES, alerts, and drop reasons
- `audit_log.json` — UTC timestamp, HitPilot/RDKit versions, config, input SHA-256, and exact counts

By default PAINS, Brenk, and NIH matches are recorded and dropped (`drop_pains`,
`drop_brenk`, and `drop_nih` all default to `true`). Pass `--no-drop-nih` (or the
equivalent flags for PAINS/Brenk) to keep those compounds and flag them only.


## Phase 2 quickstart

Phase 2 scores the compounds that passed hygiene for physicochemical
properties (including the real Ertl SA score), ChEMBL nearest-neighbor
novelty, Bemis-Murcko scaffolds, Butina clusters, and a sphere-exclusion
diverse subset.

```bash
# Fetch unique CHEMBL206 actives (pChEMBL >= 6.0, cached under data/chembl/)
.venv/bin/hitpilot fetch-chembl CHEMBL206 --pchembl-min 6.0 --max-records 2000

# Run properties, novelty, and diversity together on Phase 1 survivors
.venv/bin/hitpilot phase2 outputs/phase1_ESR1_ant/passed_compounds.csv \
  --target-chembl-id CHEMBL206 \
  --output-dir outputs/phase2_ESR1_ant

# Or run each step on its own
.venv/bin/hitpilot properties outputs/phase1_ESR1_ant/passed_compounds.csv \
  --output-dir outputs/phase2_ESR1_ant/properties
.venv/bin/hitpilot novelty outputs/phase1_ESR1_ant/passed_compounds.csv \
  --target-chembl-id CHEMBL206 --output-dir outputs/phase2_ESR1_ant/novelty
.venv/bin/hitpilot diversity outputs/phase1_ESR1_ant/passed_compounds.csv \
  --output-dir outputs/phase2_ESR1_ant/diversity
```

Unified Phase 2 outputs:

- `phase2_enriched.csv` — Phase 1 columns plus descriptors, novelty tiers, scaffolds, and cluster/sphere assignments
- `diverse_subset.csv` — sphere-exclusion representatives (no two share Morgan Tanimoto >= 0.65 by default)
- `scaffold_summary.csv` — Bemis-Murcko scaffold frequencies
- `audit_log.json` — UTC timestamp, versions, configs, input SHA-256, and exact summaries


## Phase 3 quickstart

Phase 3 adds a pluggable docking/ML rescorer CSV hook (with standard VS
enrichment metrics and a score-ordered diverse shortlist) and **ChemMCP**,
an MCP server that wraps the verified HitPilot tools.

```bash
# Leak-free baseline: train an ECFP4 random forest on LIT-PCBA _T, score _V
.venv/bin/hitpilot baseline-score ESR1_ant \
  --output outputs/phase3_ESR1_ant/baseline_ecfp4_rf_scores.csv

# Join the scorer CSV onto Phase 2 compounds, rank, and pick a diverse top-25
.venv/bin/hitpilot score outputs/phase2_ESR1_ant/phase2_enriched.csv \
  outputs/phase3_ESR1_ant/baseline_ecfp4_rf_scores.csv \
  --output-dir outputs/phase3_ESR1_ant --top-n 25

# Launch ChemMCP over stdio (or use the hitpilot-mcp console script)
.venv/bin/hitpilot mcp
# equivalent: .venv/bin/hitpilot-mcp
```

Scorer hook outputs:

- `ranked_compounds.csv` — Phase 1/2 compounds that received an external score, ranked best-first
- `shortlist_topN_diverse.csv` — score-ordered sphere-exclusion shortlist (pairwise Morgan Tanimoto < 0.65 by default)
- `audit_log.json` — SHA-256 checksums for both input files, config, enrichment metrics, and shortlist stats

Use `--lower-is-better` when the external score is a docking energy (Vina/Glide kcal/mol). The default `--higher-is-better` is for ML probabilities or pIC50 values.

Sample MCP client config:

```json
{
  "mcpServers": {
    "ChemMCP": {
      "command": ".venv/bin/hitpilot-mcp"
    }
  }
}
```


## Phase 4 quickstart

Phase 4 writes a publication-grade med-chem dossier (Markdown + a self-contained
HTML report with embedded figures) and a one-command public benchmark that
runs Phases 1-3 and then the dossier on a LIT-PCBA target.

```bash
# One-command public benchmark (reuses cached data/litpcba and data/chembl)
.venv/bin/hitpilot benchmark ESR1_ant \
  --target-chembl-id CHEMBL206 \
  --output-dir outputs/benchmark_ESR1_ant \
  --top-n 25

# Or build a dossier from existing Phase 1-3 output directories
.venv/bin/hitpilot dossier \
  --phase1-dir outputs/phase1_ESR1_ant \
  --phase2-dir outputs/phase2_ESR1_ant \
  --phase3-dir outputs/phase3_ESR1_ant \
  --output-dir outputs/dossier_ESR1_ant \
  --target ESR1_ant --target-chembl-id CHEMBL206 --top-n 25
```

Dossier outputs (under `<output_dir>/` or `outputs/benchmark_ESR1_ant/dossier/`):

- `medchem_dossier.md` — funnel, enrichment table, top-25 diverse shortlist, scaffold and alert audits, plus relative `figures/*.png` links
- `medchem_dossier.html` — the same report as a single self-contained file (figures embedded as base64 data URIs)
- `dossier_manifest.json` — SHA-256 hashes of all inputs and artifacts plus key summary metrics
- `figures/top_shortlist_grid.png`, `flagged_alerts_grid.png`, `physchem_distributions.png`, `enrichment_and_triage_funnel.png`, `scaffold_and_clusters.png`

Pass `--scores-csv` to `hitpilot benchmark` to use an external docking or ML rescorer instead of the leak-free `_T` to `_V` ECFP4 random-forest baseline.

## Five-phase roadmap

1. **Phase 1 — Hygiene and alerts.** RDKit parsing and sanitization, salt/solvent
   stripping, charge neutralization, canonical SMILES / InChIKey, exact
   deduplication, PAINS / Brenk / NIH FilterCatalog alerts, auditable logs, CLI,
   and tests.
2. **Phase 2 — Properties, novelty, and diversity.** Physicochemical filters,
   scaffold clustering, diversity selection, and novelty checks against known
   public chemistry.
3. **Phase 3 — ChemMCP, CLI, and a pluggable scorer hook.** MCP tooling plus
   an interface for external docking or ML rescoring.
4. **Phase 4 — Med-chem dossier on the public benchmark.** An end-to-end
   dossier for the chosen LIT-PCBA target, including virtual-screening
   enrichment context.
5. **Phase 5 — GitHub release and CV packaging.** Public repository, license,
   and materials ready to share.
