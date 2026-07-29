#!/usr/bin/env bash
# One-command environment setup on the GPU server.
#
#   bash scripts/setup_server.sh
#
# Creates a conda env (or falls back to a venv), installs the CPU stack and
# Boltz-2, then fetches the MMP9 receptor and builds the fragment pool.
# Safe to re-run.

set -euo pipefail
cd "$(dirname "$0")/.."
ENV_NAME="${ENV_NAME:-findaptamer}"

echo "==> 0. environment report"
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv 2>/dev/null \
  || echo "!! nvidia-smi not available — the Boltz-2 oracle needs a GPU"
python3 -V; nproc; free -g | head -2

echo
echo "==> 1. python environment"
if command -v conda >/dev/null 2>&1; then
  # shellcheck disable=SC1091
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda env list | grep -q "^${ENV_NAME} " || conda create -y -n "$ENV_NAME" python=3.11
  conda activate "$ENV_NAME"
  echo "   conda env: $ENV_NAME"
else
  [ -d .venv ] || python3 -m venv .venv
  # shellcheck disable=SC1091
  source .venv/bin/activate
  echo "   venv: .venv"
fi
python -m pip install -q --upgrade pip

echo
echo "==> 2. CPU stack"
python -m pip install -q -r requirements.txt
python -c "import RNA, numpy; print('   ViennaRNA', RNA.__version__, '| numpy', numpy.__version__)"

echo
echo "==> 3. Boltz-2 (binding oracle)"
if python -c "import boltz" 2>/dev/null; then
  echo "   boltz already installed"
else
  python -m pip install -q boltz && echo "   boltz installed" \
    || echo "   !! boltz install failed — check network/proxy, then: pip install boltz"
fi
command -v boltz >/dev/null && boltz --help >/dev/null 2>&1 \
  && echo "   boltz CLI OK" || echo "   !! boltz CLI not on PATH"

echo
echo "==> 4. MMP9 receptor (UniProt P14780 + PDB)"
python src/target/fetch_receptor.py --domain catalytic --outdir data/mmp9 \
  || echo "   !! receptor fetch failed — check outbound network, then re-run this step"

echo
echo "==> 5. fragment pool"
python src/generator/pool.py --mode generic --out data/mmp9/pool.csv

echo
echo "==> 6. self-test (CPU, no GPU needed)"
python src/target/mmp9.py >/dev/null && echo "   target module OK"
python src/oracle/interface.py
python src/scoring/benchmark_known.py >/dev/null && echo "   benchmarks OK"

cat <<'EOF'

============================================================
Setup complete. Next:

  # 1) CALIBRATE the oracle (the gate — must PASS before designing)
  python src/oracle/calibrate.py --receptor data/mmp9/mmp9_receptor.fasta \
      --out results/calibration.json

  # 2) Run the closed-loop design
  bash scripts/run_pilot.sh

============================================================
EOF
