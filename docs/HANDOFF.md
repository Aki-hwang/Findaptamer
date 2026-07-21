# HANDOFF — resume on the GPU workstation

This repo is self-contained. To continue the MMP9 aptamer project on a different
(GPU) workstation, clone it and follow this file. Nothing here depends on the dev
sandbox.

## Locked-in decisions
- **Goal: high-affinity BINDER** (diagnostic/detection), not an inhibitor.
  → No FnII activation concern; validation is MST/SPR Kd + selectivity vs MMP2/MMP7.
  → **Bar to beat: F3B ≈ 20 nM** (2′-F RNA, MMP9-selective). Aim: a DNA binder ≤ ~20 nM.
- **Method: oracle-in-the-loop** (beats AiDTA's proxy reward). Full design in
  `docs/02_improved_method.md`; critique of AiDTA in `docs/01_aidta_analysis.md`.
- **Oracle stack: Boltz-2** (MIT, commercial-OK) as the in-loop binding oracle;
  AF3 for optional consensus; MM/GBSA for the calibrated energy tier.

## What is DONE and validated (CPU, committed)
| File | Status |
|---|---|
| `src/fragments/build_library.py` | ✅ reproduces the 5,440 ss + 5,440 ds DNA fragment library |
| `src/scoring/secondary_structure.py` | ✅ DNA folding + AiDTA reward + thermodynamic metrics (ViennaRNA) |
| `src/scoring/benchmark_known.py` | ✅ folds F3B etc. → `results/benchmarks.json` |
| `src/generator/pool.py` | ✅ fragment-pool builder: generic (CPU) + HDOCK-docked-pool loader |
| `src/generator/assembler.py` | ✅ faithful AiDTA assembly engine (CPU); proves the reward accepts unstructured non-binders; WC-hairpin correctness checked over 5k assemblies |
| `src/pipeline/closed_loop.py` | ✅ generate→oracle→select driver; runs on CPU proxy, `Boltz2Oracle` drop-in; beats AiDTA-reward selection |
| `src/target/mmp9.py` | ✅ MMP9 domains, epitopes, receptors, FnII zone, benchmark aptamers |
| `src/target/receptor.py` | ✅ fetch/slice the MMP9 catalytic domain — fails loudly, never fabricates a sequence |
| `src/oracle/interface.py` | ✅ Oracle ABC + `StructureProxyOracle` (CPU; fixes AiDTA degeneracy) |
| `src/oracle/boltz2.py` | ✅ `Boltz2Oracle` GPU wrapper — ready to run, never fabricates scores |

> Note: an earlier handoff listed `assembler.py` as committed when it was not in
> the repo. It (plus `pool.py`, `pipeline/closed_loop.py`, `target/receptor.py`)
> were built and validated on CPU in the workstation session of 2026-07-21.

Quick check after clone:
```bash
pip install -r requirements.txt
python3 src/target/mmp9.py
python3 src/scoring/benchmark_known.py
python3 src/oracle/interface.py       # structured 0.80 vs unstructured(all-A) 0.39
python3 src/generator/assembler.py    # AiDTA accepts unstructured junk (oracle ~0.23)
python3 src/pipeline/closed_loop.py --rounds 12 --pop 60   # writes results/candidates.json
```

## Progress in the 2026-07-21 workstation session (CPU only — this box has no GPU)
- Installed and verified the CPU stack (ViennaRNA 2.7.2 etc.); all quick checks pass.
- Built the three missing pieces from the TODO below: `generator/pool.py`,
  `pipeline/closed_loop.py`, and `target/receptor.py`, plus the previously
  miscounted `generator/assembler.py`.
- Ran the closed loop end-to-end with `StructureProxyOracle` →
  `results/candidates.json`. Findings: (1) oracle-driven selection beats AiDTA's
  own structure-consistency selection decisively (mean top-15 ≈ 0.95 vs 0.75);
  (2) vs plain random search on the *cheap smooth CPU proxy* it is only at
  parity — expected; the sample-efficiency edge is a property of the expensive,
  rugged Boltz-2 oracle and is deliberately NOT claimed from the proxy run.
- UniProt egress (`rest.uniprot.org`) is blocked by this session's policy, so the
  real MMP9 sequence could not be fetched here and was **not** fabricated.
- **GPU check:** this session ran on a CPU-only cloud VM (4 vCPU Xeon, no GPU);
  the user's local box is a Windows machine with an **NVIDIA T400, 4 GB VRAM** —
  6× below the 24 GB minimum for the Boltz-2 MMP9+DNA co-fold (see the hardware
  table in `docs/02_improved_method.md`). **Next step is a ≥24 GB GPU** (RTX
  4090 / L40 / A5000 / A100, e.g. a cloud rental or a proper workstation), then
  `bash scripts/run_on_gpu.sh`. On Windows, use WSL2 (native `.sh`/Makefile need
  a POSIX shell) or run on a Linux GPU host.

## What to do NEXT on the workstation (in order)

> **Shortcut:** on the GPU workstation, `bash scripts/run_on_gpu.sh` does steps
> 1–5 below in one command (verifies the GPU, installs Boltz-2, fetches the MMP9
> sequence without fabricating it, smoke-tests the oracle, runs the loop →
> `results/candidates_boltz2.json`). `RECEPTOR_FASTA=P14780.fasta` if UniProt is
> unreachable. The steps below are the manual/explained version.

1. **Install GPU tools:** `pip install boltz` (downloads weights on first run).
   Optionally AF3 (non-commercial) and OpenMM/AmberTools for Stage 5.
2. **Get the MMP9 receptor sequence/structure.** On an open-network host:
   `python3 src/target/receptor.py` fetches P14780 and prints the catalytic
   domain (107-443); pass `drop_fnii=True` for the FnII-free mini construct.
   Also pull PDB **1GKC** (active-site) or **5TH6** (apo). Feed the sequence to a
   `Boltz2Config`.
3. **Smoke-test the oracle on one candidate:**
   ```python
   from src.oracle.boltz2 import Boltz2Oracle, Boltz2Config
   cfg = Boltz2Config(receptor_sequence="<MMP9 catalytic domain seq>")
   print(Boltz2Oracle(cfg).score("TCTGTCTGGGCGACATTTGCCGTGGCCACCCAGCG"))
   ```
   Confirm it returns an interface ipTM. Calibrate by also scoring F3B (should
   rank high) and a random 50-mer (should rank low).
4. **Fragment pool (Stage 1):** DONE as code. `build_generic_pool()` gives an
   oracle-targeted pool with no docking. For an AiDTA-style target-biased pool,
   HDOCK-dock the library to the epitope and load the results with
   `load_docked_pool("<hdock>.csv", higher_is_better=False)`.
5. **Close the loop (Stage 2–3) — DONE on CPU, one flag from GPU:** the driver
   `src/pipeline/closed_loop.py` already does generate → oracle → select. Flip the
   oracle:
   ```bash
   python3 src/pipeline/closed_loop.py --oracle boltz2 \
       --receptor "$(python3 src/target/receptor.py | tail -1)" --rounds 12 --pop 60
   ```
   That is the core "better than AiDTA" step, now running against real interface
   scores. (Consider raising `--rounds`, and MCTS/GFlowNet as a stronger generator.)
6. **Consensus + energy (Stage 4–5):** re-fold survivors with AF3 (multi-seed);
   short MD + MM/GBSA on the top ~100 for calibrated ΔG ranking.
7. **Shortlist ~10–30** for synthesis + MST/SPR, with selectivity vs MMP2/MMP7.

## Remaining TODO (not yet built)
- Add the FnII counter-selection `penalty_fn` (needs the Boltz-2 pose to detect
  FnII-exosite contact; the hook already exists in `ClosedLoopDesigner`).
- `src/oracle/alphafold3.py` — AF3 consensus wrapper (Stage 4).
- Stage-5 MD/MM-GBSA scripts.

## Note on git
The first two commits (151f00a, 87161b2) were authored under a personal email and
show as "Unverified" on GitHub; force-fixing them was blocked by sandbox policy.
Content is identical and intact — cosmetic only. All later commits are correctly
authored. No action needed.
