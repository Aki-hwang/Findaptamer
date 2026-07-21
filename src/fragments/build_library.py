"""
Build the nucleic-acid fragment library used as the building blocks for
docking-then-assembling aptamer design.

Reproduces the library described in Guo et al. (AiDTA, bioRxiv 2025):
  - single-stranded (ss) fragments of length 3, 4, 5, 6 nt
  - double-stranded (ds) fragments of length 3, 4, 5, 6 bp

Counts (target-independent, before any docking selection):
  ss: 4^3 + 4^4 + 4^5 + 4^6 = 5,440
  ds: 4^3 + 4^4 + 4^5 + 4^6 = 5,440   (a duplex is defined by one strand;
                                       the partner is its reverse complement)

Each fragment is emitted with its dot-bracket secondary structure in the
'&'-joined two-strand notation used by the assembler, e.g.
  ss:  ('GGT',       '...')
  ds:  ('GGC&GCC',   '(((&)))')

This is DNA (A/G/C/T). The library itself carries no target information; the
target-specific pool (44 fragments in the paper) is produced later by docking
each fragment to the protein epitope and keeping the top scorers.
"""
from __future__ import annotations
import argparse
import csv
import itertools
import json
from pathlib import Path

BASES = ("A", "G", "C", "T")
_COMP = {"A": "T", "T": "A", "G": "C", "C": "G"}


def revcomp(seq: str) -> str:
    return "".join(_COMP[b] for b in reversed(seq))


def enumerate_ss(lengths=(3, 4, 5, 6)):
    """Yield (sequence, structure) for every ss fragment of the given lengths."""
    for L in lengths:
        struct = "." * L
        for tup in itertools.product(BASES, repeat=L):
            yield "".join(tup), struct


def enumerate_ds(lengths=(3, 4, 5, 6)):
    """Yield (sequence, structure) for every fully-paired ds fragment.

    A duplex of L bp is represented as  strand1 & reverse_complement(strand1)
    with structure  ('*L) & (')'*L). Enumerating strand1 over all 4^L
    sequences enumerates every distinct Watson-Crick duplex of that length.
    """
    for L in lengths:
        struct = f"{'(' * L}&{')' * L}"
        for tup in itertools.product(BASES, repeat=L):
            s1 = "".join(tup)
            yield f"{s1}&{revcomp(s1)}", struct


def build(lengths=(3, 4, 5, 6)):
    ss = list(enumerate_ss(lengths))
    ds = list(enumerate_ds(lengths))
    return ss, ds


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", default="data/fragments", help="output directory")
    args = ap.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    ss, ds = build()

    with open(outdir / "fragments_ss.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sequence", "structure", "length", "type"])
        for seq, st in ss:
            w.writerow([seq, st, len(seq), "ss"])

    with open(outdir / "fragments_ds.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sequence", "structure", "bp", "type"])
        for seq, st in ds:
            bp = len(seq.split("&")[0])
            w.writerow([seq, st, bp, "ds"])

    summary = {
        "n_ss": len(ss),
        "n_ds": len(ds),
        "n_total": len(ss) + len(ds),
        "ss_by_len": {L: sum(1 for s, _ in ss if len(s) == L) for L in (3, 4, 5, 6)},
        "ds_by_bp": {L: sum(1 for s, _ in ds if len(s.split("&")[0]) == L) for L in (3, 4, 5, 6)},
        "mononucleotides": list(BASES),
    }
    with open(outdir / "library_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(json.dumps(summary, indent=2))
    print(f"\nWrote {outdir/'fragments_ss.csv'} and {outdir/'fragments_ds.csv'}")


if __name__ == "__main__":
    main()
