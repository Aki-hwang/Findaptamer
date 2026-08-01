#!/usr/bin/env bash
# Request a GPU allocation on the Inha AIX Slurm cluster and run the pipeline.
#
# Usage (from the login node, AFTER scripts/setup_login.sh has been run):
#
#   PARTITION=p1 EMAIL=you@inha.ac.kr bash scripts/slurm_run.sh            # interactive
#   PARTITION=p1 EMAIL=you@inha.ac.kr MODE=batch bash scripts/slurm_run.sh # sbatch
#
# Env vars:
#   PARTITION  p1 | p2 | p3        (REQUIRED — find yours with:
#                                   sacctmgr show assoc format=User,Partition where user=`whoami`)
#   EMAIL      your@inha.ac.kr     (REQUIRED — only inha.edu / inha.ac.kr accepted)
#   GPU_TYPE   a6000 | a100 | a40  (default a6000 — 48GB fits the full 337aa receptor)
#   GPU_N      1..4 (a40: 1..3)    (default 1 — scoring is serial)
#   TIME       D-HH:MM:SS          (default 1-00:00:00 = 1 day; cluster max 7-00:00:00)
#   MODE       interactive | batch (default interactive)
#   STAGE      calibrate | pose | control | pilot | both  (default both)
#              control = canonical protein-DNA positive control
#              (Zif268/1AAY) — tells you whether a failed aptamer
#              calibration means the model or the setup
#              pose = pose-reproducibility calibration (see
#              docs/06_calibration_result_ko.md); needed because the
#              scalar ipTM proved too noisy to rank candidates
#   SEEDS      predictions per sequence in calibration (default 3;
#              use 1 for a fast first check — 13 calls instead of 39)
#   NNEG       negative controls in calibration (default 12)
#
# Cluster policy notes honored here:
#   * installs/downloads happen on the login node (scripts/setup_login.sh), so the
#     GPU allocation is not held idle — the centre monitors usage and penalizes waste.
#   * the allocation is time-boxed; the pipeline checkpoints so it can resume.

set -euo pipefail
cd "$(dirname "$0")/.."

PARTITION="${PARTITION:-p2}"   # mellab -> p2 (verify: sacctmgr show assoc format=User,Partition user=$(whoami))
EMAIL="${EMAIL:-aki@inha.ac.kr}"
GPU_TYPE="${GPU_TYPE:-a6000}"
GPU_N="${GPU_N:-1}"   # the pipeline scores serially; a 2nd GPU would sit idle
TIME="${TIME:-1-00:00:00}"
MODE="${MODE:-interactive}"
STAGE="${STAGE:-both}"
SEEDS="${SEEDS:-3}"      # predictions per sequence in calibration
NNEG="${NNEG:-12}"       # negative controls in calibration
JOB="${JOB:-mmp9apt}"

case "$EMAIL" in
  *@inha.ac.kr|*@inha.edu) ;;
  *) echo "!! EMAIL must be @inha.ac.kr or @inha.edu (got: $EMAIL)"; exit 1;;
esac
if [ "$GPU_TYPE" = "a40" ] && [ "$GPU_N" -gt 3 ]; then
  echo "!! a40 allows at most 3 GPUs per request"; exit 1
fi
if [ "$GPU_N" -gt 4 ]; then echo "!! max 4 GPUs per request"; exit 1; fi

RECEPTOR="${RECEPTOR:-data/mmp9/mmp9_receptor.fasta}"

read -r -d '' PAYLOAD <<EOF || true
set -euo pipefail
# Batch stdout is a file, so Python block-buffers several KB before the
# first line ever reaches the log. That makes a running job look dead.
export PYTHONUNBUFFERED=1
cd "\$SLURM_SUBMIT_DIR"
echo "=== node: \$(hostname) ==="
nvidia-smi --query-gpu=index,name,memory.total --format=csv,noheader

# A compute node starts a fresh non-interactive shell, so conda is NOT on PATH
# even though setup_login.sh installed it — the previous check fell through to a
# .venv that does not exist and killed the job at the first line under set -e.
# Source the conda hook from its known location instead of relying on PATH.
_activated=0
for _base in "\${CONDA_HOME:-\$HOME/miniconda3}" "\$HOME/anaconda3" "\$HOME/miniforge3" "/opt/conda"; do
  if [ -f "\$_base/etc/profile.d/conda.sh" ]; then
    # shellcheck disable=SC1091
    . "\$_base/etc/profile.d/conda.sh"
    if conda activate "${ENV_NAME:-findaptamer}" 2>/dev/null; then
      _activated=1; echo "env: conda ${ENV_NAME:-findaptamer} (\$_base)"; break
    fi
  fi
done
if [ "\$_activated" = "0" ] && [ -f .venv/bin/activate ]; then
  # shellcheck disable=SC1091
  . .venv/bin/activate; _activated=1; echo "env: .venv"
fi
if [ "\$_activated" = "0" ]; then
  echo "!! could not activate the '${ENV_NAME:-findaptamer}' environment."
  echo "   Run scripts/setup_login.sh on the login node first."
  exit 1
fi
echo "python: \$(command -v python)"
command -v boltz >/dev/null || { echo "!! boltz not on PATH in this env"; exit 1; }

if ! python -c "
import sys; sys.path.insert(0,'src')
from pathlib import Path
from oracle.boltz2 import find_cached_msa
seq=''.join(l.strip() for l in Path('$RECEPTOR').read_text().splitlines() if not l.startswith('>'))
sys.exit(0 if find_cached_msa(seq) else 1)" 2>/dev/null; then
  echo; echo "=== caching receptor MSA (once per receptor; speeds up every later prediction) ==="
  python src/target/make_msa.py --receptor "$RECEPTOR" || \
    echo "   !! MSA caching failed — continuing with --use_msa_server (slower)"
fi

if [ "$STAGE" = "calibrate" ] || [ "$STAGE" = "both" ]; then
  echo; echo "=== oracle calibration (gate) ==="
  python src/oracle/calibrate.py --receptor "$RECEPTOR" \
    --seeds $SEEDS --n-negatives $NNEG --out results/calibration.json
fi
if [ "$STAGE" = "pose" ]; then
  echo; echo "=== pose-reproducibility calibration ==="
  python src/oracle/calibrate_pose.py --receptor "$RECEPTOR" \
    --seeds $SEEDS --n-negatives $NNEG --out results/calibration_pose.json
fi
if [ "$STAGE" = "control" ]; then
  echo; echo "=== canonical protein-DNA control (Zif268 / 1AAY) ==="
  python src/oracle/control_canonical.py --seeds $SEEDS \
    --out results/control_1aay.json
fi
if [ "$STAGE" = "pilot" ] || [ "$STAGE" = "both" ]; then
  echo; echo "=== closed-loop design ==="
  bash scripts/run_pilot.sh
fi
EOF

echo "==> requesting: ${GPU_N}x ${GPU_TYPE}  partition=${PARTITION}  time=${TIME}"

if [ "$MODE" = "batch" ]; then
  mkdir -p logs
  sbatch <<SBATCH
#!/usr/bin/env bash
#SBATCH --job-name=${JOB}
#SBATCH --partition=${PARTITION}
#SBATCH --gres=gpu:${GPU_TYPE}:${GPU_N}
#SBATCH --time=${TIME}
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${EMAIL}
#SBATCH --output=logs/${JOB}-%j.out
#SBATCH --error=logs/${JOB}-%j.err
#SBATCH --open-mode=append
${PAYLOAD}
SBATCH
  echo "submitted. watch with:  squeue -u \$(whoami)   /   tail -f logs/${JOB}-*.out"
else
  # interactive: the centre's documented pattern is `srun ... --pty bash`
  srun --gres=gpu:${GPU_TYPE}:${GPU_N} \
       -p "${PARTITION}" \
       --time="${TIME}" \
       -J "${JOB}" \
       --mail-type=ALL --mail-user="${EMAIL}" \
       --pty bash -lc "${PAYLOAD}; echo; echo '(allocation still open — Ctrl+D to release)'; exec bash"
fi
