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
#   GPU_N      1..4 (a40: 1..3)    (default 2)
#   TIME       D-HH:MM:SS          (default 1-00:00:00 = 1 day; cluster max 7-00:00:00)
#   MODE       interactive | batch (default interactive)
#   STAGE      calibrate | pilot | both   (default both)
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
GPU_N="${GPU_N:-2}"
TIME="${TIME:-1-00:00:00}"
MODE="${MODE:-interactive}"
STAGE="${STAGE:-both}"
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
cd "\$SLURM_SUBMIT_DIR"
echo "=== node: \$(hostname) ==="
nvidia-smi --query-gpu=index,name,memory.total --format=csv,noheader

if command -v conda >/dev/null 2>&1 && conda env list | grep -q "^${ENV_NAME:-findaptamer} "; then
  source "\$(conda info --base)/etc/profile.d/conda.sh"; conda activate "${ENV_NAME:-findaptamer}"
elif [ -d .venv ]; then source .venv/bin/activate; fi

if [ ! -s data/mmp9/mmp9_receptor_msa.a3m ] && [ ! -s data/mmp9/mmp9_receptor_msa.csv ]; then
  echo; echo "=== caching receptor MSA (once; speeds up every later prediction) ==="
  python src/target/make_msa.py --receptor "$RECEPTOR" || \
    echo "   !! MSA caching failed — continuing with --use_msa_server (slower)"
fi

if [ "$STAGE" = "calibrate" ] || [ "$STAGE" = "both" ]; then
  echo; echo "=== oracle calibration (gate) ==="
  python src/oracle/calibrate.py --receptor "$RECEPTOR" --out results/calibration.json
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
