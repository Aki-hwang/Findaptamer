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
| `src/generator/assembler.py` | ✅ faithful AiDTA assembly engine (CPU); exposes the proxy-reward degeneracy |
| `src/target/mmp9.py` | ✅ MMP9 domains, epitopes, receptors, FnII zone, benchmark aptamers |
| `src/oracle/interface.py` | ✅ Oracle ABC + `StructureProxyOracle` (CPU; fixes AiDTA degeneracy) |
| `src/oracle/boltz2.py` | ✅ `Boltz2Oracle` GPU wrapper — ready to run, never fabricates scores |

Quick check after clone:
```bash
pip install -r requirements.txt
python3 src/target/mmp9.py
python3 src/scoring/benchmark_known.py
python3 src/oracle/interface.py     # structured 0.80 vs unstructured(all-A) 0.39
```

## What to do NEXT on the workstation (in order)

1. **Install GPU tools:** `pip install boltz` (downloads weights on first run).
   Optionally AF3 (non-commercial) and OpenMM/AmberTools for Stage 5.
2. **Get the MMP9 receptor sequence/structure.** Pull UniProt **P14780** FASTA and
   PDB **1GKC** (active-site) or **5TH6** (apo, for the andecaliximab-selective
   patch). Put the catalytic-domain sequence into a `Boltz2Config`. (These were
   proxy-blocked in the sandbox; the workstation should have open network.)
3. **Smoke-test the oracle on one candidate:**
   ```python
   from src.oracle.boltz2 import Boltz2Oracle, Boltz2Config
   cfg = Boltz2Config(receptor_sequence="<MMP9 catalytic domain seq>")
   print(Boltz2Oracle(cfg).score("TCTGTCTGGGCGACATTTGCCGTGGCCACCCAGCG"))
   ```
   Confirm it returns an interface ipTM. Calibrate by also scoring F3B (should
   rank high) and a random 50-mer (should rank low).
4. **Build the MMP9 fragment pool (Stage 1):** either (a) HDOCK-dock the fragment
   library to the epitope and keep top scorers (AiDTA-style), or (b) skip docking
   and let the oracle-in-the-loop do the targeting from a generic structured pool
   (`src/generator` — a `build_generic_pool` helper is the next small TODO).
5. **Close the loop (Stage 2–3):** wire `assembler.batch_generate` (or an
   MCTS/GFlowNet generator) to `Boltz2Oracle` as the reward: generate → co-fold →
   keep top interface scores → (optionally) retrain the generator. This is the
   core "better than AiDTA" step and the main remaining build.
6. **Consensus + energy (Stage 4–5):** re-fold survivors with AF3 (multi-seed);
   short MD + MM/GBSA on the top ~100 for calibrated ΔG ranking.
7. **Shortlist ~10–30** for synthesis + MST/SPR, with selectivity vs MMP2/MMP7.

## Remaining TODO (not yet built)
- `src/generator/pool.py` — generic + HDOCK-docked fragment-pool builder.
- `src/pipeline/closed_loop.py` — generate→oracle→select driver (CPU proxy oracle
  default; Boltz2Oracle drop-in). This is the first thing to write on the workstation.
- `src/oracle/alphafold3.py` — AF3 consensus wrapper.
- Stage-5 MD/MM-GBSA scripts.

## Note on git
The first two commits (151f00a, 87161b2) were authored under a personal email and
show as "Unverified" on GitHub; force-fixing them was blocked by sandbox policy.
Content is identical and intact — cosmetic only. All later commits are correctly
authored. No action needed.
