#!/usr/bin/env bash
# Findaptamer — one-command GPU runbook (Stage 3: Boltz-2 in-loop binding oracle).
#
# Run THIS on the actual GPU workstation (not the Claude Code web sandbox, which
# is CPU-only). It installs deps, obtains the MMP9 receptor sequence WITHOUT
# fabricating it, smoke-tests the Boltz-2 oracle, then runs the closed loop with
# the real binding oracle and writes results/candidates_boltz2.json.
#
# Usage:
#   bash scripts/run_on_gpu.sh                 # fetch MMP9 seq from UniProt
#   RECEPTOR_FASTA=P14780.fasta bash scripts/run_on_gpu.sh   # offline: use a local FASTA
#   ROUNDS=16 POP=80 bash scripts/run_on_gpu.sh              # override loop size
#
# Requirements on the workstation: NVIDIA GPU (>=24 GB for Boltz-2 DNA co-fold),
# Python 3.10+, ~50 GB free disk (Boltz weights + MSA cache), open outbound
# network (Boltz weights download + ColabFold MSA server; UniProt for the seq).
set -euo pipefail
cd "$(dirname "$0")/.."
ROUNDS="${ROUNDS:-12}"; POP="${POP:-60}"; TOP="${TOP:-20}"

echo "==> 0. Verify a GPU is actually present"
if ! command -v nvidia-smi >/dev/null 2>&1; then
  echo "ERROR: nvidia-smi not found. This host has no NVIDIA GPU / driver."
  echo "       Boltz-2 needs a GPU. Run this on the GPU workstation." >&2
  exit 1
fi
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || { echo "ERROR: GPU not usable" >&2; exit 1; }

echo "==> 1. Install the stack (CPU deps + Boltz-2)"
python3 -m pip install -q -r requirements.txt
python3 -m pip install -q boltz    # first predict() downloads weights to ~/.boltz

echo "==> 2. Obtain the MMP9 catalytic-domain sequence (never fabricated)"
if [[ -n "${RECEPTOR_FASTA:-}" && -f "${RECEPTOR_FASTA}" ]]; then
  RECEPTOR="$(python3 -c "import sys; sys.path.insert(0,'src'); \
from target.receptor import catalytic_domain_from_fasta as f; \
print(f(open('${RECEPTOR_FASTA}').read()))")"
  echo "    using local FASTA ${RECEPTOR_FASTA} (${#RECEPTOR} aa)"
else
  # receptor.py fetches P14780 and prints the catalytic domain on the last line;
  # if egress is blocked it exits non-zero with instructions (no guessed seq).
  if ! RECEPTOR="$(python3 src/target/receptor.py | tail -1)"; then
    echo "ERROR: could not fetch the MMP9 sequence and none was provided." >&2
    echo "       Download it and re-run with RECEPTOR_FASTA=P14780.fasta:" >&2
    echo "       curl -o P14780.fasta https://rest.uniprot.org/uniprotkb/P14780.fasta" >&2
    exit 1
  fi
  echo "    fetched P14780 catalytic domain (${#RECEPTOR} aa)"
fi
case "$RECEPTOR" in *[!ACDEFGHIKLMNPQRSTVWY]*|"") echo "ERROR: bad receptor sequence" >&2; exit 1;; esac

echo "==> 3. Smoke-test the Boltz-2 oracle (calibration: F3B-like should outrank a random 50-mer)"
python3 - "$RECEPTOR" <<'PY'
import sys; sys.path.insert(0, "src")
from oracle.boltz2 import Boltz2Oracle, Boltz2Config
rec = sys.argv[1]
o = Boltz2Oracle(Boltz2Config(receptor_sequence=rec))
import random
probe = {
  "candidate":  "TCTGTCTGGGCGACATTTGCCGTGGCCACCCAGCG",
  "random50":   "".join(random.Random(0).choice("AGCT") for _ in range(50)),
}
for name, seq in probe.items():
    s = o.score(seq)
    print(f"    {name:10s} interface_score={s.value:.4f}  {s.details}")
PY

echo "==> 4. Run the closed loop with the REAL binding oracle"
python3 src/pipeline/closed_loop.py --oracle boltz2 --receptor "$RECEPTOR" \
    --rounds "$ROUNDS" --pop "$POP" --top "$TOP" --out results
# rename so the CPU-proxy candidates.json is not overwritten
if [[ -f results/candidates.json ]]; then mv -f results/candidates.json results/candidates_boltz2.json; fi

echo "==> DONE. Binding-driven shortlist -> results/candidates_boltz2.json"
echo "    Next (Stage 4-5): AF3 multi-seed consensus + MD/MM-GBSA on the top ~100,"
echo "    then synthesize ~10-30 and validate by MST/SPR (selectivity vs MMP2/MMP7)."
