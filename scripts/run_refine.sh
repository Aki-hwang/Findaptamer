#!/usr/bin/env bash
# Stages 4-5: consensus validation + MM/GBSA energy ranking + final report.
#
# Run AFTER scripts/run_pilot.sh has produced results/pilot/top_candidates.csv.
#
#   bash scripts/run_refine.sh
#   TOP=20 SEEDS=3 USE_CHAI=1 NS=1.0 bash scripts/run_refine.sh
#
# Env:
#   TOP       candidates carried into consensus      (default 20)
#   SEEDS     Boltz-2 seeds per candidate            (default 3)
#   USE_CHAI  1 = add Chai-1 as a second model       (default 0)
#   ENERGY_N  candidates carried into MM/GBSA        (default 10)
#   NS        MD nanoseconds per candidate           (default 0.5)
#   EPITOPE   report overlap with a named epitope    (default: none)

set -euo pipefail
cd "$(dirname "$0")/.."

RECEPTOR="${RECEPTOR:-data/mmp9/mmp9_receptor.fasta}"
CANDS="${CANDS:-results/pilot/top_candidates.csv}"
TOP="${TOP:-20}"
SEEDS="${SEEDS:-3}"
USE_CHAI="${USE_CHAI:-0}"
ENERGY_N="${ENERGY_N:-10}"
NS="${NS:-0.5}"
EPITOPE="${EPITOPE:-}"

if command -v conda >/dev/null 2>&1 && conda env list | grep -q "^${ENV_NAME:-findaptamer} "; then
  # shellcheck disable=SC1091
  source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate "${ENV_NAME:-findaptamer}"
elif [ -d .venv ]; then
  # shellcheck disable=SC1091
  source .venv/bin/activate
fi

echo "==> preflight"
[ -f "$CANDS" ] || { echo "!! missing $CANDS — run scripts/run_pilot.sh first"; exit 1; }
[ -f "$RECEPTOR" ] || { echo "!! missing $RECEPTOR"; exit 1; }
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || { echo "!! no GPU"; exit 1; }

echo
echo "==> Stage 4: consensus (top $TOP, $SEEDS seeds$([ "$USE_CHAI" = 1 ] && echo ' + Chai-1'))"
python src/pipeline/consensus.py \
  --receptor "$RECEPTOR" \
  --candidates "$CANDS" \
  --top "$TOP" \
  --seeds "$SEEDS" \
  ${USE_CHAI:+$([ "$USE_CHAI" = 1 ] && echo --use-chai)} \
  ${EPITOPE:+--epitope "$EPITOPE"} \
  --out results/consensus

echo
echo "==> Stage 5: MM/GBSA energy (top $ENERGY_N, ${NS} ns each)"
if python -c "import openmm" 2>/dev/null; then
  python src/pipeline/energy.py \
    --consensus results/consensus/consensus.json \
    --top "$ENERGY_N" --ns "$NS" --only-pass --out results/energy \
    || echo "   !! energy stage failed — continuing to the report with consensus only"
else
  echo "   openmm not installed — skipping. To enable: pip install openmm"
fi

echo
echo "==> Final report"
python src/pipeline/report.py --n "$TOP" --out results/report

cat <<'EOF'

============================================================
Done.
  results/consensus/consensus_ranked.csv   consensus-validated candidates
  results/energy/energy_ranked.csv         MM/GBSA dG ranking
  results/report/REPORT.md                 synthesis shortlist (read this)
============================================================
EOF
