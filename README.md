# Findaptamer — de novo ML design of DNA aptamers against MMP9

Goal: computationally design DNA aptamers that bind (and ideally inhibit) human
**MMP9** (matrix metalloproteinase-9, gelatinase B; UniProt **P14780**), using a
machine-learning pipeline that is **more accurate than the AiDTA proxy-reward
approach** (Guo et al., bioRxiv 2025) by putting a real structure + affinity
oracle *in the loop*.

## Why not just re-run AiDTA

AiDTA's reinforcement-learning reward is secondary-structure self-consistency —
it never sees the protein and never optimizes binding; binding is only filtered
post-hoc and open-loop. See `docs/01_aidta_analysis.md` for the full critique.
Our thesis: close that loop with a modern protein–nucleic-acid co-folding /
affinity model. See `docs/02_improved_method.md`.

## Repository layout

```
src/
  fragments/build_library.py     # exact 5,440 ss + 5,440 ds DNA fragment library
  scoring/secondary_structure.py # DNA folding + reward (ViennaRNA, DNA params)
  target/mmp9.py                 # MMP9 domains, epitopes, benchmarks
  target/receptor.py             # fetch/slice MMP9 catalytic domain (never fabricates)
  generator/pool.py              # fragment-pool builder (generic + HDOCK-docked loader)
  generator/assembler.py         # faithful AiDTA assembly engine (exposes its degeneracy)
  oracle/interface.py            # Oracle ABC + CPU StructureProxyOracle
  oracle/boltz2.py               # GPU Boltz-2 binding oracle (drop-in)
  pipeline/closed_loop.py        # generate -> oracle -> select driver (the core step)
docs/                            # method analysis & design
data/fragments/                  # library summary
results/                         # ranked candidates + reports
```

## Resuming on a GPU workstation

**Read `docs/HANDOFF.md` first.** It lists the locked-in decisions (goal =
high-affinity binder; bar to beat = F3B ~20 nM), what is already validated on CPU,
and the exact next steps (install Boltz-2, get the MMP9 receptor, close the loop).

## Status

CPU-validated (committed; all run without a GPU):
- [x] Fragment library reproduces the paper exactly (10,880 fragments).
- [x] DNA secondary-structure scoring via ViennaRNA (DNA Mathews 2004).
- [x] Critical analysis of AiDTA locating the improvement target.
- [x] MMP9 target/epitope module (P14780; receptors, epitopes, FnII zone, benchmarks).
- [x] Oracle interface + CPU `StructureProxyOracle` (fixes AiDTA degeneracy).
- [x] `Boltz2Oracle` GPU wrapper (ready to run; never fabricates scores).
- [x] Fragment-pool builder — generic (CPU) + HDOCK-docked-pool loader (`generator/pool.py`).
- [x] Faithful AiDTA assembly engine — proves AiDTA's reward accepts unstructured
      non-binders (`generator/assembler.py`); WC-hairpin correctness checked over 5k assemblies.
- [x] **Closed-loop generate→oracle→select driver** (`pipeline/closed_loop.py`) — runs
      end-to-end on CPU with the proxy oracle; `Boltz2Oracle` is a one-line drop-in.
      Beats AiDTA-reward selection decisively (see `results/candidates.json`).
- [x] Receptor helper (`target/receptor.py`) — fetches/slices the MMP9 catalytic
      domain when the network allows; fails loudly, never fabricates a sequence.

Needs the GPU workstation (see `docs/HANDOFF.md`):
- [ ] Swap `Boltz2Oracle` into the closed loop (turns the CPU-validated machinery
      into a binding-driven design run). This is now a `--oracle boltz2` flag flip.
- [ ] Real MMP9 receptor sequence (UniProt egress was policy-blocked in this
      session; `target/receptor.py` pulls it on an open-network host).
- [ ] Optional: HDOCK-dock the fragment library to the epitope for a target-biased pool.
- [ ] AF3 consensus + MD/MM-GBSA energy tier (Stages 4–5).
- [ ] Shortlist → synthesis → MST/SPR validation.

## Environment notes

Runs so far are CPU-only (this sandbox has no GPU). The heavy stages
(co-folding / affinity oracle, physics refinement) need a GPU workstation; the
method doc specifies exactly which stage needs what.
