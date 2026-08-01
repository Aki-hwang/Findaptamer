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
RESUME="${RESUME:-1}"          # 1 = continue from checkpoint; 0 = fresh start
MAX_HOURS="${MAX_HOURS:-0}"    # stop+checkpoint before the Slurm --time expires
CALIB="${CALIB:-results/calibration.json}"

# --- environment activation -------------------------------------------------
# Compute nodes start a fresh non-interactive shell where conda is not on PATH,
# so `command -v conda` fails and a naive .venv fallback aborts the run under
# `set -e`. Source the conda hook from its install location instead.
activate_env() {
  local name="${ENV_NAME:-findaptamer}" base
  for base in "${CONDA_HOME:-$HOME/miniconda3}" "$HOME/anaconda3" "$HOME/miniforge3" "/opt/conda"; do
    if [ -f "$base/etc/profile.d/conda.sh" ]; then
      # shellcheck disable=SC1091
      . "$base/etc/profile.d/conda.sh"
      if conda activate "$name" 2>/dev/null; then
        echo "env: conda $name ($base)"; return 0
      fi
    fi
  done
  if [ -f .venv/bin/activate ]; then
    # shellcheck disable=SC1091
    . .venv/bin/activate; echo "env: .venv"; return 0
  fi
  echo "!! could not activate '$name' — run scripts/setup_login.sh first"; return 1
}
activate_env || exit 1

echo "==> preflight"
[ -f "$RECEPTOR" ] || { echo "!! missing $RECEPTOR — run: python src/target/fetch_receptor.py"; exit 1; }
[ -f "$POOL" ]     || { echo "!! missing $POOL — run: python src/generator/pool.py"; exit 1; }
command -v boltz >/dev/null || { echo "!! boltz not found — run scripts/setup_login.sh"; exit 1; }
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || { echo "!! no GPU"; exit 1; }

if [ -f "$CALIB" ]; then
  if python - "$CALIB" <<'PY'
import json, sys
s = json.load(open(sys.argv[1]))["summary"]
au, z = s.get("auroc", 0), s.get("separation_z")
ok = (s.get("oracle") == "boltz2_interface" and au >= 0.75
      and (z is None or z >= 1.0))
print(f"   calibration: oracle={s.get('oracle')} AUROC={au} z={z} "
      f"(positive sequences={s.get('n_positive_sequences')}) -> {'PASS' if ok else 'FAIL'}")
if s.get("evidence_note"):
    print(f"   note: {s['evidence_note']}")
sys.exit(0 if ok else 1)
PY
  then :; else
    if [ "${FORCE:-0}" = "1" ]; then
      echo "!! Calibration did not pass, but FORCE=1 was set - proceeding anyway."
      echo "   The oracle has not been shown to separate binders from negatives;"
      echo "   treat every score from this run as unvalidated."
    else
      echo "!! Calibration did not pass. Fix the oracle before designing."
      echo "   Re-run:  python src/oracle/calibrate.py --receptor $RECEPTOR"
      echo "   Options: try the other receptor slice (--domain catalytic_nofn),"
      echo "            raise --seeds, or inspect results/calibration.json."
      echo "   To proceed anyway (results will be unvalidated): FORCE=1 $0"
      exit 1
    fi
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
  --max-hours "$MAX_HOURS" \
  $([ "$RESUME" = "1" ] && echo --resume) \
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
