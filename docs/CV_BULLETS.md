# CV and portfolio language for HitPilot

Use these as-is or lightly adapt them. Every number is taken from the verified
HitPilot 0.1.0 run on public LIT-PCBA ESR1_ant / ChEMBL CHEMBL206 data. Do not
round them into a stronger claim, and do not add unpublished models or private
datasets.

## CV bullets

- Built **HitPilot**, an open-source, auditable hit-triage and medicinal-chemistry copilot (Python, RDKit, scikit-learn, MCP): structure hygiene, official PAINS/Brenk/NIH alerts, physicochemical filters, ChEMBL novelty, scaffold diversity, a pluggable docking/ML rescorer hook, and a shareable med-chem dossier, released under MIT at [github.com/hentablihamza/hitpilot](https://github.com/hentablihamza/hitpilot).
- Delivered a fully reproducible end-to-end triage of the public **LIT-PCBA ESR1_ant** library (PubChem AID 743080): **5,050** input compounds standardized and screened down to **1,696** hygiene-passed structures (40 actives / 1,656 inactives) with machine-readable SHA-256 audit logs.
- Implemented leak-aware virtual-screening evaluation on official LIT-PCBA train/validation splits, catching and removing **405** standardized SMILES overlaps between `_T` and `_V` before scoring so enrichment is not inflated by train–test structure leakage.
- Engineered a pluggable CSV hook for *any* docking or ML rescorer (IDs or SMILES; higher-is-better or docking energies) that reports sklearn ROC-AUC/PR-AUC, standard EF at 1/5/10%, and official RDKit BEDROC, then builds a score-ordered diverse shortlist (Morgan Tanimoto &lt; 0.65). A leak-free ECFP4 random-forest baseline on the hygiene-matched held-out set (N = 368, 6 actives) gave **ROC-AUC 0.580**, **EF1% 15.33**, and **BEDROC(α=20) 0.209**.
- Computed nearest-neighbor novelty against **1,101** unique ChEMBL **CHEMBL206** actives (pChEMBL ≥ 6) and mapped the passed library into **822** unique Bemis–Murcko scaffolds and **1,275** Butina clusters, with a top-25 diverse shortlist for med-chem review.
- Shipped **ChemMCP**, a 9-tool MCP server so any MCP-compatible agent can run verified RDKit standardization, descriptors, ChEMBL/PubChem lookups, diversity selection, rescorer evaluation, and dossier generation over stdio.
- Packaged the work for public reuse: **65** automated tests, Ruff-clean code, GitHub Actions CI on Python 3.11 and 3.12, and a publication-grade dossier (Markdown + self-contained HTML with embedded structure grids and enrichment plots).

## LinkedIn / portfolio showcase

I built and open-sourced HitPilot, a reproducible hit-triage and med-chem copilot for public virtual-screening data. On the LIT-PCBA ESR1_ant benchmark it takes 5,050 public compounds through RDKit hygiene (1,696 passed), ChEMBL CHEMBL206 novelty (1,101 reference actives; 822 Murcko scaffolds; 1,275 Butina clusters), and a pluggable docking/ML rescorer hook that leak-checks the official train/validation split (405 overlapping structures removed) before reporting ROC-AUC, enrichment factors, and BEDROC and writing a shareable dossier. A bundled ECFP4 random-forest baseline on the hygiene-matched held-out set reached ROC-AUC 0.580 and EF1% 15.33 — not a claim of a state-of-the-art model, but a transparent, auditable demonstration that any external score CSV can be evaluated the same way. ChemMCP exposes nine verified chemistry tools to MCP-compatible agents. Code, tests (65), and example figures: https://github.com/hentablihamza/hitpilot

## Elevator pitch

HitPilot is an open-source copilot that turns a public screening library into an auditable med-chem shortlist and dossier — without fabricating chemistry. It standardizes molecules, applies official structural alerts, checks novelty against ChEMBL, measures diversity, and lets you drop in any docking or ML score file to get leak-checked enrichment and a ready-to-review report. Computational chemists, med chemists, and LLM agents can all use the same verified RDKit tools.
