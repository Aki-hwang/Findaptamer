#!/usr/bin/env bash
# Setup on the LOGIN NODE (no GPU allocation held).
#
#   bash scripts/setup_login.sh
#
# The centre monitors GPU usage and penalizes waste, so everything that does not
# need a GPU — env creation, pip installs, Boltz weight download — happens here.
# Afterwards, request a GPU with scripts/slurm_run.sh and go straight to compute.

set -euo pipefail
cd "$(dirname "$0")/.."
ENV_NAME="${ENV_NAME:-findaptamer}"

echo "==> 0. where am I"
hostname
echo "   (this should be the login node, e.g. login-n1 — do NOT run heavy compute here)"

echo
echo "==> 1. my Slurm partition"
sacctmgr show assoc format=User,Partition where user=`whoami` 2>/dev/null \
  || echo "   !! sacctmgr unavailable — ask the centre which partition (p1/p2/p3) you are in"

echo
echo "==> 2. python environment"
# HPC nodes often lack python3-venv and users have no sudo, so fall back to a
# user-space Miniconda install (self-contained, no root required).
CONDA_HOME="${CONDA_HOME:-$HOME/miniconda3}"
if ! command -v conda >/dev/null 2>&1 && [ -x "$CONDA_HOME/bin/conda" ]; then
  # shellcheck disable=SC1091
  source "$CONDA_HOME/etc/profile.d/conda.sh"
fi

if command -v conda >/dev/null 2>&1; then
  # shellcheck disable=SC1091
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda env list | grep -q "^${ENV_NAME} " || conda create -y -n "$ENV_NAME" python=3.11
  conda activate "$ENV_NAME"
  echo "   conda env: $ENV_NAME ($(python -V))"
elif python3 -m venv --help >/dev/null 2>&1 && python3 -c "import ensurepip" 2>/dev/null; then
  [ -d .venv ] || python3 -m venv .venv
  # shellcheck disable=SC1091
  source .venv/bin/activate
  echo "   venv: .venv ($(python -V))"
else
  echo "   no conda and python3-venv is unavailable (no sudo on a shared node)."
  echo "   installing Miniconda into $CONDA_HOME (user-space, no root needed)..."
  MC=/tmp/miniconda_$USER.sh
  curl -fsSL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o "$MC" \
    || wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O "$MC"
  bash "$MC" -b -p "$CONDA_HOME"
  rm -f "$MC"
  # shellcheck disable=SC1091
  source "$CONDA_HOME/etc/profile.d/conda.sh"
  "$CONDA_HOME/bin/conda" init bash >/dev/null 2>&1 || true
  conda create -y -n "$ENV_NAME" python=3.11
  conda activate "$ENV_NAME"
  echo "   Miniconda installed; conda env: $ENV_NAME ($(python -V))"
fi
python -m pip install -q --upgrade pip

echo
echo "==> 3. CPU stack"
python -m pip install -q -r requirements.txt
python -c "import RNA, numpy; print('   ViennaRNA', RNA.__version__, '| numpy', numpy.__version__)"

echo
echo "==> 4. Boltz-2"
python -c "import boltz" 2>/dev/null && echo "   already installed" \
  || python -m pip install -q boltz && echo "   installed"

echo
echo "==> 5. pre-download Boltz weights (do it HERE, not on a GPU allocation)"
python - <<'PY' || echo "   !! weight pre-download failed; it will download on first GPU run"
try:
    from boltz.main import download_boltz2, MODEL_URL  # noqa: F401
except Exception:
    pass
from pathlib import Path
cache = Path.home() / ".boltz"
cache.mkdir(exist_ok=True)
print(f"   boltz cache: {cache}")
try:
    from boltz.main import download_boltz2 as dl
    dl(cache); print("   weights ready")
except Exception as e:
    print(f"   (could not pre-fetch programmatically: {e})")
    print("   fallback: run one tiny prediction on the login node is NOT allowed;")
    print("   weights will download automatically on the first GPU job instead.")
PY

echo
echo "==> 6. receptor + fragment pool (vendored sequence, no network needed)"
python src/target/fetch_receptor.py --domain catalytic --outdir data/mmp9 --pdb
python src/generator/pool.py --mode generic --out data/mmp9/pool.csv

echo
echo "==> 7. CPU self-test"
python src/oracle/interface.py
python src/scoring/benchmark_known.py >/dev/null && echo "   benchmarks OK"

cat <<'EOF'

============================================================
Login-node setup complete. Now request a GPU and run:

  PARTITION=<p1|p2|p3> EMAIL=<you@inha.ac.kr> bash scripts/slurm_run.sh

  # or non-interactive:
  PARTITION=p1 EMAIL=you@inha.ac.kr MODE=batch bash scripts/slurm_run.sh

Defaults: 2x a6000, 1 day. Override with GPU_TYPE / GPU_N / TIME.
============================================================
EOF
