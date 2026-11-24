# Changelog

All notable changes to HitPilot are documented in this file.

## 0.1.0 — 2026-11-24

First public release. The original five-phase roadmap is complete, and all
chemistry is computed on public LIT-PCBA, ChEMBL, and PubChem data with
verified RDKit and scikit-learn code.

### Phase 1 — Hygiene and alerts

- RDKit parsing and sanitization, salt/solvent stripping, charge neutralization,
  canonical SMILES / InChIKey, and exact structure deduplication.
- Official RDKit FilterCatalog alerts for PAINS (A/B/C), Brenk, and NIH, with
  drop-or-flag-only switches and a machine-readable audit log (SHA-256 of the
  input, versions, and exact drop counts).
- LIT-PCBA fetcher and loader, plus a `hitpilot hygiene` CLI.

### Phase 2 — Properties, novelty, and diversity

- RDKit descriptors, the official Ertl SA score, and Lipinski / Veber flags.
- Morgan ECFP4 nearest-neighbor novelty against cached or live ChEMBL target
  actives (default CHEMBL206), with novelty tiers and optional PubChem PUG REST
  lookups.
- Bemis–Murcko scaffolds, Butina clustering, and score-aware sphere-exclusion
  diversity selection (`hitpilot phase2` and the step-wise subcommands).

### Phase 3 — ChemMCP and pluggable scorer hook

- CSV hook for any docking or ML rescorer (`higher_is_better` or docking
  energies), with ID join, canonical SMILES fallback, unmatched-ID reporting,
  and a score-ordered diverse shortlist.
- Standard VS enrichment: sklearn ROC-AUC / PR-AUC, EF at 1/5/10%, and RDKit
  `CalcBEDROC` at α = 20 and α = 80.5.
- Leak-free LIT-PCBA baseline: an ECFP4 random forest trained only on `_T` and
  scored only on `_V` after dropping standardized SMILES shared with train
  (405 overlaps removed on ESR1_ant).
- ChemMCP stdio server exposing verifiable hit-triage tools, plus
  `hitpilot baseline-score`, `hitpilot score`, and `hitpilot mcp`.

### Phase 4 — Med-chem dossier and public benchmark

- Publication-grade figures (structure grids, physchem plots, triage funnel,
  ROC and enrichment, scaffolds and clusters).
- Markdown dossier and a self-contained HTML report with figures embedded as
  base64 data URIs, plus `dossier_manifest.json` hashes.
- One-command `hitpilot benchmark ESR1_ant` that runs Phases 1–4 and reuses
  cached public data.

### Phase 5 — GitHub release and CV packaging

- MIT license, contributing guide, CITATION.cff, changelog, and a rewritten
  README with architecture diagram, worked ESR1_ant example figures, and
  audience-specific guidance.
- GitHub Actions CI (Ruff + pytest) on Python 3.11 and 3.12, issue/PR
  templates, and committed example assets under `docs/images/` and `examples/`.

### Verified public benchmark (LIT-PCBA ESR1_ant / CHEMBL206)

- 5,050 input compounds; 3,906 unique standardized structures; 1,696
  hygiene-passed (40 actives / 1,656 inactives).
- 1,101 unique ChEMBL CHEMBL206 reference actives; 822 unique Murcko scaffolds;
  1,275 Butina clusters.
- 405 `_T`/`_V` standardized SMILES overlaps removed before scoring.
- Leak-free, hygiene-matched held-out set (N = 368, 6 actives): ROC-AUC 0.580,
  EF 1% 15.33, BEDROC (α = 20) 0.209 for the bundled ECFP4 random-forest
  baseline. These figures demonstrate the hook; they are not a claim of a
  state-of-the-art rescorer.
