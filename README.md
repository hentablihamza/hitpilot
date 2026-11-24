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
