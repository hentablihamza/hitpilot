# Contributing to HitPilot

Thanks for helping improve an auditable, public-data hit-triage toolbox.
HitPilot is deliberately conservative: chemistry must be computed with real
RDKit code or live public APIs (LIT-PCBA, ChEMBL, PubChem). Never invent
SMILES, scores, IDs, or enrichment metrics, and never add unpublished models
or private datasets to the default pipeline.

## Development setup

Python 3.11+ is required. The project is developed with [uv](https://docs.astral.sh/uv/):

```bash
# Clone and enter the repo
git clone https://github.com/hentablihamza/hitpilot.git
cd hitpilot

# Create a virtualenv and install the package plus dev tools
uv venv --python 3.12
uv pip install --python .venv/bin/python -e ".[dev]"

# Always run tools through the project virtualenv
.venv/bin/pytest -v
.venv/bin/ruff check src tests
.venv/bin/hitpilot --help
```

`pip install -e ".[dev]"` works the same way if you already have a virtualenv.
Raw LIT-PCBA and ChEMBL downloads live under `data/` and run outputs under
`outputs/`; both directories are gitignored. Cached public files from a prior
run are reused automatically by `hitpilot benchmark`.

## Code style

- Line length 100, linted with Ruff (`E`, `F`, `I`, `W`) as configured in `pyproject.toml`.
- Public functions and Pydantic models should stay typed and documented.
- Prefer extending an existing Phase 1–4 module over introducing a parallel path.
- Tests go in `tests/` and must run offline by default (mock HTTP). Mark any
  live-network check with `@pytest.mark.network` so GitHub Actions can skip it.

## Adding a new structural alert catalog

Alert screening lives in `src/hitpilot/hygiene.py` and uses RDKit
`FilterCatalog`. To add a catalog (for example another official RDKit set):

1. Extend `HygieneConfig` in `src/hitpilot/models.py` with an `enable_*` flag
   and a corresponding `drop_*` flag (default drop behavior should be explicit).
2. Register the catalog in `build_filter_catalog()` using the official
   `FilterCatalogParams.FilterCatalogs` member — do not invent SMARTS.
3. Map matches onto a stable drop reason and catalog name in `collect_alerts()`.
4. Add unit tests in `tests/test_hygiene.py` that use real molecules known to
   trigger (or not trigger) the official catalog, plus a flag-only CLI test if
   you expose `--no-drop-*`.
5. Document the new flag in the README command reference.

## Adding a scorer hook (docking or ML CSV)

The pluggable hook is `src/hitpilot/scoring.py`. An external rescorer is just a
CSV with a compound ID (or SMILES) and a numeric score:

1. Follow `ScorerHookConfig`: `id_column`, optional `smiles_column` fallback,
   `score_column`, and `higher_is_better` (true for probabilities / pIC50, false
   for docking energies such as Vina or Glide kcal/mol).
2. Keep ranking, merge-count reporting, and leak-free evaluation in
   `merge_external_scores()` and `compute_enrichment_metrics()`. Enrichment
   must use `sklearn.metrics` for ROC-AUC / PR-AUC, the standard VS definition
   of EF, and `rdkit.ML.Scoring.Scoring.CalcBEDROC` for BEDROC.
3. Diverse shortlists must reuse the Phase 2 Morgan fingerprint + sphere
   exclusion engine so pairwise Tanimoto stays below the cutoff.
4. Add coverage in `tests/test_scoring.py` for both score directions, unmatched
   IDs, SMILES fallback joins, and (if you change metrics) exact numeric
   checks on a synthetic ranking.
5. Wire any new CLI options through `hitpilot score` / `hitpilot benchmark`
   (`--scores-csv` already accepts an arbitrary rescorer table).

## Adding a ChemMCP tool

The MCP server is `src/hitpilot/mcp_server.py` (`mcp` 2.x `MCPServer`).

1. Implement the chemistry in a testable package function first (hygiene,
   properties, novelty, diversity, scoring, or dossier). The MCP tool should
   be a thin, typed wrapper that returns JSON-serializable dicts.
2. Register the tool with `@server.tool()` inside `create_chem_mcp_server()`.
3. Add the tool name to `EXPECTED_TOOLS` in `tests/test_mcp_server.py` and
   cover the input schema plus at least one in-process invocation (mock any
   live ChEMBL or PubChem call).
4. Document the tool in the README ChemMCP table and, if it is a CLI-facing
   pipeline, in the command reference.
5. Reinstall the console script after entry-point changes:

   ```bash
   uv pip install --python .venv/bin/python -e . --no-deps
   ```

## Pull requests

Please include tests, keep `.venv/bin/ruff check src tests` clean, and say
which public data or official RDKit behavior you used. Do not commit
`data/`, `outputs/`, or `.venv/`.
