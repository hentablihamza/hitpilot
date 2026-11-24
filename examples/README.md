# Examples

This directory holds a one-command demo and a committed snapshot of the
HitPilot 0.1.0 med-chem dossier for the public LIT-PCBA `ESR1_ant`
benchmark (ChEMBL target CHEMBL206).

## One-command demo

From the repository root, with the project virtualenv installed as in the
[README](../README.md) or [CONTRIBUTING](../CONTRIBUTING.md) guide:

```bash
bash examples/run_demo.sh
```

The script runs:

```bash
hitpilot benchmark ESR1_ant \
  --target-chembl-id CHEMBL206 \
  --output-dir outputs/benchmark_ESR1_ant \
  --top-n 25
```

The first run downloads public LIT-PCBA and ChEMBL files into `data/` (gitignored)
and writes Phase 1–4 outputs under `outputs/benchmark_ESR1_ant/` (also gitignored).
Repeat runs reuse the cached public data. Open
`outputs/benchmark_ESR1_ant/dossier/medchem_dossier.html` in a browser for the
self-contained report.

To score your own docking or ML table instead of the bundled leak-free ECFP4
baseline, pass `--scores-csv path/to/scores.csv` (see `hitpilot score --help`).

## Committed snapshot

- [`example_dossier_ESR1_ant.md`](example_dossier_ESR1_ant.md) — Markdown dossier
  from the verified 0.1.0 run. Three figures are stored in
  [`docs/images/`](../docs/images/) so GitHub visitors can view results without
  running the pipeline. The remaining two figures (alert examples and scaffold
  plot) are produced by the demo command above.
- SHA-256 paths inside the snapshot refer to files a local benchmark run writes
  under `outputs/`; those raw tables are not committed.

Verified headline metrics from that run are tabulated in the README.
