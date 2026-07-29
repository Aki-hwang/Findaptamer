#!/usr/bin/env bash
# One-command pilot design run on the GPU server.
#
#   bash scripts/run_pilot.sh                 # default: 8 rounds x 96 candidates
#   ROUNDS=12 PER_ROUND=128 bash scripts/run_pilot.sh
#   OUT=results/run2 bash scripts/run_pilot.sh
#
# Refuses to start unless the oracle calibration has PASSED — driving generation
# with an oracle that cannot separate known binders from negatives would just
# optimize noise.

set -euo pipefail
cd "$(dirname "$0")/.."

RECEPTOR="${RECEPTOR:-data/mmp9/mmp9_receptor.fasta}"
POOL="${POOL:-data/mmp9/pool.csv}"
ROUNDS="${ROUNDS:-8}"
PER_ROUND="${PER_ROUND:-96}"
MIN_LEN="${MIN_LEN:-30}"
MAX_LEN="${MAX_LEN:-45}"
TOP="${TOP:-30}"
OUT="${OUT:-results/pilot}"
CALIB="${CALIB:-results/calibration.json}"

# activate the env created by setup_server.sh
if command -v conda >/dev/null 2>&1 && conda env list | grep -q "^${ENV_NAME:-findaptamer} "; then
  # shellcheck disable=SC1091
  source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate "${ENV_NAME:-findaptamer}"
elif [ -d .venv ]; then
  # shellcheck disable=SC1091
  source .venv/bin/activate
fi

echo "==> preflight"
[ -f "$RECEPTOR" ] || { echo "!! missing $RECEPTOR — run: python src/target/fetch_receptor.py"; exit 1; }
[ -f "$POOL" ]     || { echo "!! missing $POOL — run: python src/generator/pool.py"; exit 1; }
command -v boltz >/dev/null || { echo "!! boltz not found — run scripts/setup_server.sh"; exit 1; }
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || { echo "!! no GPU"; exit 1; }

if [ -f "$CALIB" ]; then
  if python - "$CALIB" <<'PY'
import json, sys
s = json.load(open(sys.argv[1]))["summary"]
ok = s.get("auroc", 0) >= 0.75 and s.get("oracle") == "boltz2_interface"
print(f"   calibration: oracle={s.get('oracle')} AUROC={s.get('auroc')} -> {'PASS' if ok else 'FAIL'}")
sys.exit(0 if ok else 1)
PY
  then :; else
    echo "!! Calibration did not pass. Fix the oracle before designing."
    echo "   Re-run: python src/oracle/calibrate.py --receptor $RECEPTOR"
    exit 1
  fi
else
  echo "!! No calibration found at $CALIB."
  echo "   Run first: python src/oracle/calibrate.py --receptor $RECEPTOR"
  exit 1
fi

echo
echo "==> closed-loop design  (rounds=$ROUNDS x $PER_ROUND, len $MIN_LEN-$MAX_LEN)"
python src/pipeline/closed_loop.py \
  --oracle boltz2 \
  --receptor "$RECEPTOR" \
  --pool "$POOL" \
  --rounds "$ROUNDS" \
  --per-round "$PER_ROUND" \
  --min-len "$MIN_LEN" \
  --max-len "$MAX_LEN" \
  --prefilter \
  --top "$TOP" \
  --out "$OUT"

echo
echo "============================================================"
echo "Done. Results:"
echo "  $OUT/top_candidates.csv   ranked candidate aptamers"
echo "  $OUT/result.json          full record incl. per-round history"
echo
echo "Next (Stage 4-5): consensus re-fold of the top candidates with a second"
echo "model + multiple seeds, then short MD + MM/GBSA for calibrated dG ranking."
echo "============================================================"
