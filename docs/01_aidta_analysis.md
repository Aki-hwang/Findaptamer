# AiDTA — method summary and critical analysis

Source: Guo, Guo, Qian, He, Qian, Wang, Huang. *De novo design of protein-binding
aptamers through deep reinforcement learning assembly of nucleic acid fragments.*
bioRxiv 2025.06.01.657174 (preprint). Code: github.com/Fudan-HQLab/AiDTA.

## 1. What AiDTA does (pipeline)

1. **Fragment library (target-independent).** All single-stranded DNA fragments of
   length 3–6 nt (5,440) and all fully-paired double-stranded fragments of 3–6 bp
   (5,440). 3D structures predicted with 3dRNA/DNA, energy-minimized with
   ParmBSC1 (AmberTools23).
2. **Fragment docking (target-specific).** Define the target epitope from
   antibody–antigen complex structures (residues within 6 Å). Globally dock every
   fragment to the protein with **HDOCK**. Keep fragments whose best pose centers
   within 6 Å of the epitope; rank by docking score; take top 5 ss and top 5 ds per
   length → 40, plus the 4 mononucleotides → a **44-fragment pool**.
3. **Assembly by RL (AiDTA agent).** AlphaZero-style: state = partially assembled
   sequence encoded as a 6×L matrix; action = (choose a fragment) × (choose a
   linkage site); a policy–value CNN (conv + 7 residual blocks + policy/value
   heads) guides **MCTS** (100 simulations per step, PUCT). Assemble until length
   ≥ L_min = 50 nt.
4. **Reward.** Fold the assembled sequence (RNAstructure `Fold`), compute the
   per-nucleotide paired/unpaired agreement between the folded structure and the
   structure implied by the constituent fragments. Similarity ≥ 0.9 → reward +1,
   else −1. ~500 epochs → ~20,000 aptamers with similarity ≥ 0.9.
5. **Post-hoc structural filtering (the actual binding selection).**
   3dRNA/DNA 3D models → geometric filter (5′–3′ Euclidean & geodesic span must
   exceed the epitope span) → ~15,000 → **HDOCK** global docking, nucleotide–epitope
   distance < 7 Å → ~8,000 → top 100 by HDOCK re-docked with **HADDOCK** → top 20 →
   **AlphaFold3** complex prediction → ~6 confirmed at the epitope → synthesize.
6. **Experimental validation.** MST (Cy5-labeled aptamers, Monolith X). Reported
   Kd in the nanomolar range; best per target 84.9 / 93.0–94.4 / 226 / 375 nM
   across BA.2.86 RBD, B.1.1.529 RBD, TL1A, vWF-A1.

## 2. What is genuinely good about it

- **True de novo**: needs no SELEX data and no predetermined aptamer structure —
  more general than RaptGen / AptaDiff / RhoDesign.
- **Site-specific**: design is focused on a chosen epitope.
- **Sample-efficient at the wet-lab end**: < 10 synthesized per target.
- The RL assembler produces diverse, low-homology, well-structured sequences
  (Shannon entropy ≈ 1.8/2.0, mean pairwise Levenshtein > 28).

## 3. The core weakness (where we can beat it)

**The generative model never optimizes binding.** Its RL reward is *secondary-
structure self-consistency* — "does the assembled aptamer fold into a structure
that agrees with the shapes of the fragments it was built from." That is a proxy
for "the fragments keep their docked binding modes," and it is only loosely
coupled to actual affinity. Concretely:

- **The reward contains no protein.** Two sequences with identical fold-consistency
  are indistinguishable to the agent even if one binds MMP9 tightly and the other
  not at all. The agent learns to make *foldable* sequences, not *binders*.
- **Binding selection is entirely post-hoc and open-loop.** 20,000 → 6 is done by
  docking/AF3 filters that never feed back into the generator. The generator cannot
  learn from a near-miss; it just keeps sampling and the filter discards.
- **Fragment scoring is the weakest kind of docking.** Rigid global docking of
  3–6-mers is noisy and promiscuous; tiny fragments have shallow, non-specific
  interfaces and HDOCK's statistical potential was not built for them. Pool quality
  (and thus everything downstream) rests on this.
- **Cascade of surrogates, each with error, none calibrated to affinity.**
  HDOCK score → HADDOCK score → AF3 ipTM are three different rankings; the paper
  keeps whatever survives all three, but none is an affinity estimate, so the final
  ranking that reaches the bench is only weakly predictive of Kd. The hit rate
  (≈ 3–6 of ~7 tested < 1 µM) reflects this.
- **No affinity model anywhere.** Nothing in the pipeline predicts Kd; ranking is by
  geometric/interface heuristics.

## 4. Design implications for an improved MMP9 pipeline

The improvement thesis is simple: **put a real structure+affinity oracle in the
loop.** Replace (or augment) the fold-consistency proxy reward with a modern
protein–nucleic-acid co-folding model that (a) predicts the aptamer–MMP9 complex
and (b) scores interface confidence / affinity, and **feed that signal back into
generation** (active learning / RL / iterative refinement) instead of using it only
as a terminal filter. This directly optimizes the objective AiDTA only approximates.

Foundations already built in this repo and reused unchanged:
- `src/fragments/build_library.py` — the exact 5,440 + 5,440 fragment library.
- `src/scoring/secondary_structure.py` — DNA folding via ViennaRNA (DNA Mathews
  2004), reproducing AiDTA's fold-consistency reward **and** adding thermodynamic
  robustness metrics (ensemble frequency, diversity, positional entropy) the proxy
  reward ignores.

The oracle, the target definition, and the exact loop are finalized in
`docs/02_improved_method.md` once the target and tool survey complete.
