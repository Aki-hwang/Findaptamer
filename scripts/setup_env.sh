#!/usr/bin/env bash
# Findaptamer — turnkey environment setup.
#
# Creates an isolated Python venv (.venv), installs the CPU stack, and — when a
# GPU is present (or FORCE_BOLTZ=1) — installs Boltz-2. Then runs the CPU
# regression checks so you KNOW the environment works before spending GPU time.
#
# Run once on any machine (workstation or laptop):
#   bash scripts/setup_env.sh
#   source .venv/bin/activate
# Then:  bash scripts/run_on_gpu.sh        (GPU box)     — real binding oracle
#   or:  make run-cpu                       (any box)     — CPU proxy loop
set -euo pipefail
cd "$(dirname "$0")/.."
PYTHON="${PYTHON:-python3}"

echo "==> Python: $($PYTHON --version 2>&1)"
"$PYTHON" -c 'import sys; assert sys.version_info>=(3,10), "need Python >=3.10"'

echo "==> Creating venv at .venv"
"$PYTHON" -m venv .venv
# shellcheck disable=SC1091
source .venv/bin/activate
python -m pip install -q --upgrade pip

echo "==> Installing CPU stack (requirements.txt)"
python -m pip install -q -r requirements.txt

GPU=0
if command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi >/dev/null 2>&1; then GPU=1; fi
if [[ "$GPU" == "1" || "${FORCE_BOLTZ:-0}" == "1" ]]; then
  echo "==> GPU detected (or FORCE_BOLTZ=1): installing Boltz-2"
  python -m pip install -q boltz
else
  echo "==> No GPU here: skipping Boltz-2 (huge CUDA build, useless on CPU)."
  echo "    The CPU proxy loop still runs. Install Boltz-2 on the GPU box:"
  echo "    FORCE_BOLTZ=1 bash scripts/setup_env.sh   (or: pip install boltz)"
fi

echo "==> Verifying the environment (CPU regression checks)"
python src/fragments/build_library.py           >/dev/null && echo "    [ok] fragment library"
python src/scoring/benchmark_known.py            >/dev/null && echo "    [ok] benchmark folding"
python src/oracle/interface.py                   >/dev/null && echo "    [ok] proxy oracle"
python src/generator/assembler.py                >/dev/null && echo "    [ok] assembler"
python src/pipeline/closed_loop.py --rounds 4 --pop 24 --out /tmp/findaptamer_verify >/dev/null \
  && echo "    [ok] closed loop"    # --out keeps the curated results/candidates.json intact

echo ""
echo "==> Environment ready. GPU present: $([[ $GPU == 1 ]] && echo yes || echo NO)"
echo "    activate:  source .venv/bin/activate"
if [[ "$GPU" == "1" ]]; then
  echo "    run real:  bash scripts/run_on_gpu.sh"
else
  echo "    run CPU:   make run-cpu    (Boltz-2 needs a GPU box)"
fi
