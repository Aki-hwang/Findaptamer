"""
Characterize known MMP9 aptamers as design baselines.

Folds each benchmark aptamer with the appropriate ViennaRNA parameter set
(RNA vs DNA) and reports MFE, structure, GC, and ensemble decisiveness. This
establishes what "a real MMP9 binder" looks like structurally so generated
candidates can be compared to prior art rather than scored in a vacuum.
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

import RNA

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from target.mmp9 import BENCHMARK_APTAMERS  # noqa: E402


def _fold(seq: str, kind: str) -> dict:
    md = RNA.md()
    if "DNA" in kind.upper():
        RNA.params_load_DNA_Mathews2004()
    else:
        RNA.params_load_RNA_Turner2004()
    fc = RNA.fold_compound(seq, md)
    ss, mfe = fc.mfe()
    fc.pf()
    gc = (seq.count("G") + seq.count("C")) / len(seq)
    return {
        "length": len(seq),
        "mfe": round(float(mfe), 2),
        "structure": ss,
        "gc_fraction": round(gc, 3),
        "ens_freq_mfe": round(float(fc.pr_structure(ss)), 4),
    }


def main():
    out = []
    for apt in BENCHMARK_APTAMERS:
        row = {k: apt[k] for k in ("name", "type", "affinity_nM", "function")}
        if apt.get("sequence"):
            row["fold"] = _fold(apt["sequence"], apt["type"])
        else:
            row["fold"] = "sequence not publicly available"
        out.append(row)
    Path("results").mkdir(exist_ok=True)
    Path("results/benchmarks.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
