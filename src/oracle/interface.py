"""
Oracle interface for the closed-loop generator.

An Oracle scores a candidate DNA aptamer for its fitness as an MMP9 binder.
The whole point of this project (vs. AiDTA) is that the *generator is driven by
an oracle*, so the oracle is a swappable component:

  * StructureProxyOracle  — CPU-only. A structural-fitness PRIOR (not a binding
                            predictor). Unlike AiDTA's reward it is NOT gameable
                            by unstructured sequences: it explicitly rewards a
                            decisive, stable, appropriately-paired fold and
                            penalizes the all-unpaired degenerate optimum. Use it
                            to pre-filter and to run the loop where no GPU exists.

  * Boltz2Oracle          — GPU. The real binding oracle: co-folds aptamer+MMP9
                            and scores the interface (ipTM / interface-PAE) at the
                            target epitope. This is the accuracy edge over AiDTA.
                            (see src/oracle/boltz2.py)

Both return an OracleScore with a single `.value` in [0,1] (higher = better) plus
a `details` dict, so the generator/pipeline code is oracle-agnostic.
"""
from __future__ import annotations
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scoring.secondary_structure import fold_stability, ensemble_metrics  # noqa: E402


@dataclass
class OracleScore:
    value: float                     # scalar fitness in [0,1], higher is better
    details: dict = field(default_factory=dict)
    kind: str = ""


class Oracle(ABC):
    @abstractmethod
    def score(self, sequence: str) -> OracleScore:
        ...

    def score_batch(self, sequences):
        return [self.score(s) for s in sequences]


def _clamp(x, lo=0.0, hi=1.0):
    return max(lo, min(hi, x))


class StructureProxyOracle(Oracle):
    """CPU structural-fitness prior.

    Combines, for a DNA aptamer:
      * fold decisiveness   : ensemble frequency of the MFE structure (want high)
      * fold stability      : MFE per nucleotide vs a reference scale (want stable)
      * structuredness      : paired fraction in a target band (want ~0.35-0.75;
                              penalizes BOTH unstructured and fully-duplexed)
      * GC sanity           : GC fraction in a synthesizable/foldable band
    The structuredness term is what removes AiDTA's degenerate all-unpaired optimum.
    Weights are exposed for tuning/active learning.
    """

    def __init__(self, weights=None, paired_band=(0.35, 0.75), gc_band=(0.35, 0.70),
                 mfe_ref=-0.35):
        self.w = weights or {"decisive": 0.35, "stable": 0.25,
                             "structured": 0.30, "gc": 0.10}
        self.paired_band = paired_band
        self.gc_band = gc_band
        self.mfe_ref = mfe_ref   # kcal/mol per nt considered "well stabilized"

    @staticmethod
    def _band_score(x, lo, hi):
        """1.0 inside [lo,hi], decaying linearly outside over one band-width."""
        if lo <= x <= hi:
            return 1.0
        width = (hi - lo) or 1e-6
        d = (lo - x) if x < lo else (x - hi)
        return _clamp(1.0 - d / width)

    def score(self, sequence: str) -> OracleScore:
        seq = sequence.replace("&", "")
        st = fold_stability(seq)
        en = ensemble_metrics(seq)
        decisive = _clamp(en["ens_freq_mfe"])                 # already ~[0,1]
        stable = _clamp(st["mfe_per_nt"] / self.mfe_ref)      # >=1 if very stable
        stable = min(stable, 1.0)
        structured = self._band_score(st["paired_fraction"], *self.paired_band)
        gc = self._band_score(st["gc_fraction"], *self.gc_band)
        w = self.w
        value = (w["decisive"] * decisive + w["stable"] * stable +
                 w["structured"] * structured + w["gc"] * gc)
        return OracleScore(
            value=round(value, 4),
            kind="structure_proxy",
            details={
                "decisive": round(decisive, 3), "stable": round(stable, 3),
                "structured": round(structured, 3), "gc": round(gc, 3),
                "mfe": st["mfe"], "paired_fraction": round(st["paired_fraction"], 3),
                "gc_fraction": round(st["gc_fraction"], 3),
                "ens_freq_mfe": round(en["ens_freq_mfe"], 3),
                "mfe_structure": st["mfe_structure"],
            },
        )


if __name__ == "__main__":
    o = StructureProxyOracle()
    # F3B-like structured vs an all-A unstructured control (AiDTA degenerate case)
    for name, s in [("structured", "TCTGTCTGGGCGACATTTGCCGTGGCCACCCAGCGAGAGCTGGGTGGCCAC"),
                    ("unstructured_AAAA", "A" * 50)]:
        sc = o.score(s)
        print(f"{name:20s} value={sc.value}  structured={sc.details['structured']}"
              f"  paired={sc.details['paired_fraction']}  mfe={sc.details['mfe']}")
