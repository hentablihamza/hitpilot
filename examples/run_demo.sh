#!/usr/bin/env bash
# One-command public benchmark demo for HitPilot (LIT-PCBA ESR1_ant / CHEMBL206).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT}"

if [[ -x "${ROOT}/.venv/bin/hitpilot" ]]; then
  HITPILOT="${ROOT}/.venv/bin/hitpilot"
elif command -v hitpilot >/dev/null 2>&1; then
  HITPILOT="hitpilot"
else
  echo "hitpilot is not installed. From the repo root run:" >&2
  echo "  uv venv --python 3.12" >&2
  echo "  uv pip install --python .venv/bin/python -e \".[dev]\"" >&2
  exit 1
fi

exec "${HITPILOT}" benchmark ESR1_ant \
  --target-chembl-id CHEMBL206 \
  --output-dir outputs/benchmark_ESR1_ant \
  --top-n 25
