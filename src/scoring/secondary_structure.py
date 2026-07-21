"""
DNA secondary-structure scoring for aptamer candidates.

This module is the local, fully-reproducible substrate that AiDTA obtains from
RNAstructure's `Fold`/`ct2dot`. We use ViennaRNA with the DNA
(Mathews 2004) energy parameters so that folding is appropriate for
single-stranded DNA aptamers rather than RNA.

It provides:
  * fold_dna(seq)                 -> (mfe_structure, mfe_energy)
  * structure_similarity(a, b)    -> paired/unpaired agreement in [0,1]
                                     (identical definition to AiDTA's reward:
                                      per-position match of paired '(' ')' vs
                                      unpaired '.', averaged over length)
  * ensemble_metrics(seq)         -> partition-function based robustness:
                                     ensemble MFE frequency, ensemble diversity,
                                     mean positional entropy
  * fold_stability(seq)           -> MFE, normalized MFE (per nt), GC fraction
  * score_candidate(seq, ...)     -> a single dict of all metrics

`structure_similarity` lets us reproduce the AiDTA assembly reward: compare the
folded structure of an assembled aptamer to the "assembled" structure implied
by its constituent fragments. A value >= 0.9 is a successful assembly (reward +1
in AiDTA). We keep that convention available but also expose the richer,
thermodynamically-grounded metrics that the paper's proxy reward ignores.
"""
from __future__ import annotations
from functools import lru_cache
import RNA

_DNA_LOADED = False


def _ensure_dna_params():
    global _DNA_LOADED
    if not _DNA_LOADED:
        RNA.params_load_DNA_Mathews2004()
        _DNA_LOADED = True


def _strip(seq: str) -> str:
    return seq.replace("&", "")


@lru_cache(maxsize=200_000)
def fold_dna(seq: str):
    """Return (mfe_structure, mfe_energy_kcal_per_mol) for a DNA sequence."""
    _ensure_dna_params()
    s = _strip(seq)
    fc = RNA.fold_compound(s, RNA.md())
    ss, mfe = fc.mfe()
    return ss, float(mfe)


def _binary(struct: str) -> str:
    """Map dot-bracket to paired(1)/unpaired(0), ignoring '&'."""
    return (
        struct.replace("&", "")
        .replace(".", "0")
        .replace("(", "1")
        .replace(")", "1")
    )


def structure_similarity(struct_a: str, struct_b: str) -> float:
    """Per-position paired/unpaired agreement in [0,1] (AiDTA reward definition).

    Both structures must correspond to the same-length sequence. Returns the
    fraction of positions where both are paired or both are unpaired.
    """
    a, b = _binary(struct_a), _binary(struct_b)
    n = min(len(a), len(b))
    if n == 0:
        return 0.0
    return sum(1 for i in range(n) if a[i] == b[i]) / n


def assembly_reward(seq: str, assembled_structure: str, threshold: float = 0.9):
    """AiDTA-style reward: fold the sequence, compare to the structure implied by
    the constituent fragments, return (+1/-1, similarity, folded_structure)."""
    folded, _ = fold_dna(seq)
    sim = structure_similarity(folded, assembled_structure)
    return (1 if sim >= threshold else -1), sim, folded


def ensemble_metrics(seq: str) -> dict:
    """Partition-function metrics describing how robust the fold is.

    ens_freq_mfe : Boltzmann probability of the MFE structure (higher = the
                   aptamer folds decisively into one shape).
    ens_diversity: mean base-pair distance across the ensemble (lower = less
                   conformational heterogeneity).
    mean_entropy : mean positional (Shannon) entropy of the pairing (lower =
                   more well-defined structure).
    """
    _ensure_dna_params()
    s = _strip(seq)
    fc = RNA.fold_compound(s, RNA.md())
    ss, _ = fc.mfe()
    fc.pf()
    ens_freq = float(fc.pr_structure(ss))
    try:
        ens_div = float(fc.mean_bp_distance())
    except Exception:
        ens_div = float("nan")
    # positional entropy from base-pair probability matrix
    try:
        bpp = fc.bpp()  # 1-indexed matrix
        n = len(s)
        import math
        entropies = []
        for i in range(1, n + 1):
            p_paired = 0.0
            for j in range(1, n + 1):
                if i < j:
                    p_paired += bpp[i][j]
                elif j < i:
                    p_paired += bpp[j][i]
            p_un = max(0.0, 1.0 - p_paired)
            h = 0.0
            for p in (p_paired, p_un):
                if p > 1e-9:
                    h -= p * math.log2(p)
            entropies.append(h)
        mean_entropy = sum(entropies) / n if n else float("nan")
    except Exception:
        mean_entropy = float("nan")
    return {
        "ens_freq_mfe": ens_freq,
        "ens_diversity": ens_div,
        "mean_entropy": mean_entropy,
    }


def fold_stability(seq: str) -> dict:
    s = _strip(seq)
    ss, mfe = fold_dna(seq)
    gc = (s.count("G") + s.count("C")) / len(s) if s else 0.0
    paired = ss.count("(")
    return {
        "length": len(s),
        "mfe": mfe,
        "mfe_per_nt": mfe / len(s) if s else float("nan"),
        "gc_fraction": gc,
        "paired_fraction": (2 * paired / len(s)) if s else 0.0,
        "mfe_structure": ss,
    }


def score_candidate(seq: str, assembled_structure: str | None = None) -> dict:
    """Full local score for one aptamer candidate."""
    out = fold_stability(seq)
    out.update(ensemble_metrics(seq))
    if assembled_structure is not None:
        r, sim, folded = assembly_reward(seq, assembled_structure)
        out["assembly_similarity"] = sim
        out["assembly_reward"] = r
    return out


if __name__ == "__main__":
    demo = "TCTGTCTGGGCGACATTTGCCGTGGCCACCCAGCGAGAGCTGGGTGGCCAC"
    import json
    print(json.dumps(score_candidate(demo), indent=2, default=str))
