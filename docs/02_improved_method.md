# Improved MMP9 aptamer-design pipeline — "oracle-in-the-loop"

Goal: design de novo **DNA** aptamers against **MMP9** that are **more accurate
than AiDTA** by closing the loop AiDTA leaves open — replacing its
secondary-structure-similarity proxy reward with a real protein–nucleic-acid
**structure + energy oracle**, and adding cross-model consensus and a calibrated
physics tier.

This design is grounded in a 2024–2026 tool survey (see references at end). Two
findings shape everything:

1. **No validated one-shot ML Kd predictor exists for aptamer–protein pairs**
   (Boltz-2's affinity head is small-molecule-only; AF3 ipTM is a *weak*
   aptamer-screening discriminator per ACS Synth. Biol. 2025). So accuracy comes
   from a **cascade + consensus**, not one magic score.
2. **MMP9-specific trap:** nucleic acids that bind the FnII exosite *activate*
   MMP9 (Shimada 2018). The pipeline must **counter-select FnII** and design
   toward the andecaliximab allosteric patch (inhibitor) or the S1′ cleft.

## The bar to beat

| Reference binder | Type | Kd | Function |
|---|---|---|---|
| AiDTA best (other targets) | DNA | ~85–95 nM | binder |
| **F3B** | 2′-F RNA | **~20 nM** | MMP9 binder (imaging), selective vs MMP2/7 |
| **A3** | ssDNA 40-mer | **~13 nM** | ⚠ *activator* — do NOT reproduce this mode |
| LVMH G4 series | ssDNA G-quadruplex | (patent) | **inhibitor** (closest prior art) |

Success = a **DNA** aptamer to the **inhibitory** epitope with Kd competitive
with ~13–20 nM **and** confirmed inhibition (not activation).

## Pipeline

```
Stage 0  Target prep + epitope         CPU        MMP9 catalytic domain, +Zn/Ca,
         definition                                protein MSA (once), hotspots
Stage 1  Fragment library + epitope     CPU        10,880 ss/ds DNA fragments →
         docking (AiDTA front end)                 HDOCK to epitope → pool
Stage 2  Generator (assembly)           small GPU  MCTS/AlphaZero *or* GFlowNet
         with a REAL reward                        over the fragment action space
Stage 3  Fast in-loop oracle            A100/H100  Boltz-2 (or Chai-1) co-fold
         (the core fix)                            aptamer+MMP9; reward = interface
                                                    ipTM + interface-PAE at hotspots,
                                                    minus FnII-contact penalty
Stage 4  Consensus co-fold (top ~2k)    A100/H100  AF3 + Boltz-2 (+Chai-1), multi-seed;
                                          80 GB     keep pose-consistent, confident hits
Stage 5  Physics energy tier            RTX4090/   short explicit-solvent MD (OpenMM/
         (top ~100)                      A100      Amber) + MM/GBSA/PBSA ΔG ranking
Stage 6  Active learning (optional)     CPU/GPU    synthesize 10–30, MST/SPR Kd +
                                          + lab     inhibition assay → retrain oracle
Stage 7  Validation                     wet lab    MST/SPR Kd; zymography/FRET for
                                                    INHIBITION (rule out activation)
```

### Why each stage beats AiDTA
- **Stage 2–3 replace the proxy reward with binding.** AiDTA rewards
  fold-consistency (which, as our `assembler.py` demo shows, is *gameable* — an
  unstructured sequence scores 1.0). Here the reward is a co-folded interface
  score *at the MMP9 epitope*, so the generator learns to make **binders**, and
  the FnII penalty steers away from the activation trap AiDTA has no concept of.
- **Stage 4 adds consensus.** AiDTA trusts a single HDOCK→HADDOCK→AF3 cascade;
  documented failure modes (HADDOCK3 aptamer anomalies; AF3 weak aptamer
  screening) mean single scores mislead. Multi-model/multi-seed agreement is the
  cheapest antidote.
- **Stage 5 adds a calibrated energy.** AiDTA computes **no** binding energy;
  MM/GBSA/PBSA is the best-validated quantitative ranker for nucleic-acid
  complexes and converts "confident pose" into "ranked ΔG".
- **Stage 6 closes the loop with ground truth** — the only proven route past
  current oracle noise.

## Hardware requirements (answers "what workstation?")

| Need | Minimum | Comfortable | Used by |
|---|---|---|---|
| **In-loop oracle (Boltz-2, DNA co-fold, pose-only)** | 1× GPU **24 GB** (RTX 4090 / L40 / A5000) | 1× **A100 40 GB** | Stage 3 (dominant cost) |
| **Consensus AF3 for nucleic-acid complexes** | 1× **A100/H100 80 GB** (AF3 needs 80 GB for NA) | H100 80 GB | Stage 4 |
| **MD + MM/GBSA** | RTX 4090 24 GB | A100 | Stage 5 |
| CPU stages (MSA, fragment docking, MM/GBSA post) | 16–32 cores | 64 cores | Stage 0,1,5 |
| System RAM | 64 GB | 128 GB | all |
| Disk | 500 GB | 1–2 TB (AF3 DBs ~630 GB) | Stage 0/4 |

**Minimal viable workstation:** 1× RTX 4090 (24 GB) + 64 GB RAM + 1 TB disk runs
Stages 0–3 and 5 (Boltz-2 pose-only + MD/MM-GBSA), i.e. the whole improved loop
using **Boltz-2 as the sole oracle** (MIT license, commercial-safe).
**Ideal:** add **1× A100/H100 80 GB** to run AF3 consensus (Stage 4) and larger
Boltz-2 batches. This sandbox has **no GPU**, so Stages 3–5 cannot run here.

## Licensing note
- **Boltz-1/Boltz-2, Chai-1, Protenix** — open weights, **commercial-OK** (MIT /
  Apache-2.0). Preferred for a translational project.
- **AlphaFold3** — weights on request, **non-commercial only**. Use for
  consensus/benchmarking; do not build a commercial deliverable solely on it.

## What runs in THIS sandbox now (CPU-only, no GPU)
- Stage 0 target/epitope definition — done (`src/target/mmp9.py`); receptor
  sequence fetch/slice — done (`src/target/receptor.py`, blocked egress handled).
- Fragment library — done (`src/fragments/build_library.py`).
- DNA folding / structure scoring — done (`src/scoring/secondary_structure.py`).
- Fragment pool (generic + docked loader) — done (`src/generator/pool.py`).
- Assembly generator — done (`src/generator/assembler.py`); its demo proves the
  AiDTA reward accepts unstructured non-binders (the degeneracy this method fixes).
- **Closed loop generate→oracle→select** — done (`src/pipeline/closed_loop.py`),
  run end-to-end with the CPU `StructureProxyOracle`. It beats AiDTA-reward
  selection decisively; against random search on the cheap smooth proxy it is at
  parity (expected — the sample-efficiency win belongs to the expensive Boltz-2
  oracle and is not claimed from the proxy). Output: `results/candidates.json`.
- Boltz-2 / AF3 oracle wrappers — authored as ready-to-run interfaces
  (`src/oracle/`); `Boltz2Oracle` is a `--oracle boltz2` drop-in for the loop,
  executed once a GPU workstation is attached.

## References
AiDTA bioRxiv 2025.06.01.657174 · AF3 Nature 2024 (s41586-024-07487-w), weights
non-commercial · Boltz-2 bioRxiv 2025.06.14.659707 (MIT), affinity = small-mol
only · Chai-1 (Apache-2.0) · AF3 weak aptamer screening: ACS Synth. Biol. 2025
(5c00196) · HADDOCK3 aptamer failure modes: bioRxiv 2026.05.11.724398 · MM/PBSA
aptamer ΔG: Sci. Rep. 2025 (s41598-025-12186-1) · MMP9 FnII potentiation: Biochem.
J. 2018 475:1597 · andecaliximab epitope: J Biol Chem 2017 PMID 28235803 (5TH6/5TH9)
· F3B: PLOS ONE 2016 (pone.0149387) · LVMH G4 inhibitors: US9902961B2.
