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
  target/                        # MMP9 receptor & epitope definition
  generator/                     # aptamer generator (assembly / RL)
  oracle/                        # in-the-loop structure+affinity oracle
docs/                            # method analysis & design
data/mmp9/                       # target structures & sequences
results/                         # ranked candidates + reports
configs/                         # run configs
```

## Status

Foundations (target-independent, validated in this environment):
- [x] Fragment library reproduces the paper exactly (10,880 fragments).
- [x] DNA secondary-structure scoring via ViennaRNA (DNA Mathews 2004) — validated
      against AiDTA's example sequence.
- [x] Critical analysis of AiDTA locating the improvement target.
- [ ] MMP9 target/epitope module (research in progress).
- [ ] In-the-loop oracle + generator (compute-dependent; see method doc).

## Environment notes

Runs so far are CPU-only (this sandbox has no GPU). The heavy stages
(co-folding / affinity oracle, physics refinement) need a GPU workstation; the
method doc specifies exactly which stage needs what.
