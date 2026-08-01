"""
G-quadruplex analysis of the published MMP-9 aptamers.

Motivation. Both binding oracles we built failed (docs/06): the co-folding
score was noise-dominated and the predicted pose was irreproducible, while the
identical pipeline scored a canonical protein-DNA complex at ipTM 0.974 +/-
0.0007. So we cannot compute "does this sequence bind MMP-9". What we CAN do is
use a prior extracted from experiments that already happened.

Three independent SELEX campaigns against MMP-9 -- an RNA one (8F14A, G-quartet
confirmed by CD/Tm), a DNA biosensor one (MMP9-DNA-30), and a patent DNA one
(LVMH, "G-quadruplex structure" claimed) -- converged on G-quadruplex folds.
This module quantifies that convergence instead of asserting it, using two
published, deterministic methods:

  * G4Hunter (Bedrat, Lacroix & Mergny, Nucleic Acids Res 2016;44:1746) --
    per-base G-richness/G-skewness, averaged over a window. |score| >= 1.2 is
    the authors' threshold for likely G4 formation.
  * The G4 loop motif regex, in both the 3-tetrad form (G3+ columns) and the
    2-tetrad form (G2+ columns) -- the latter is required because TBA-class
    aptamers form only two tetrads and are invisible to the 3-tetrad pattern.

And one control that has to be run before any of this is used for design:
local alignment against the two most famous G4 aptamers. G4 aptamers recur
across unrelated SELEX campaigns, so a "novel MMP-9 G4" that is really the
thrombin aptamer with a new 5' tag would be a selection artefact, not a hit.
Better to find that out here than after synthesis.

Usage:
    python src/analysis/g4.py                 # analyse the curated set
    python src/analysis/g4.py --seq ACGT...   # analyse arbitrary sequences
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from target.mmp9_aptamers import (  # noqa: E402
    KNOWN_APTAMERS, REFERENCE_G4_APTAMERS,
)

# G4Hunter's own threshold for "likely forms a G4" (NAR 2016, Fig. 5).
G4HUNTER_THRESHOLD = 1.2

# Loop lengths 1-7 nt is the standard permissive window used by QGRS/G4Hunter
# companion tools. Two variants because tetrad count differs by aptamer class.
MOTIF_3TETRAD = re.compile(r"(?=(G{3,}\w{1,7}G{3,}\w{1,7}G{3,}\w{1,7}G{3,}))")
MOTIF_2TETRAD = re.compile(r"(?=(G{2,}\w{1,7}G{2,}\w{1,7}G{2,}\w{1,7}G{2,}))")


def _runs(seq: str):
    """Yield (base, start, length) for each homopolymer run."""
    i = 0
    while i < len(seq):
        j = i
        while j < len(seq) and seq[j] == seq[i]:
            j += 1
        yield seq[i], i, j - i
        i = j


def g4hunter_profile(seq: str) -> list[int]:
    """Per-base G4Hunter score.

    A run of k consecutive G scores +min(k, 4) on every base of the run; a run
    of k consecutive C scores -min(k, 4); A/T/U score 0. Capping at 4 is what
    keeps a long G-tract from dominating -- four is all a tetrad column needs.
    """
    seq = seq.upper().replace("U", "T")
    prof = [0] * len(seq)
    for base, start, k in _runs(seq):
        if base == "G":
            v = min(k, 4)
        elif base == "C":
            v = -min(k, 4)
        else:
            continue
        for p in range(start, start + k):
            prof[p] = v
    return prof


def g4hunter_score(seq: str) -> float:
    prof = g4hunter_profile(seq)
    return sum(prof) / len(prof) if prof else 0.0


def g4hunter_windows(seq: str, window: int = 25):
    """Max windowed score — a G4 in one half of a long sequence is diluted by
    the whole-sequence mean, so the windowed maximum is the honest readout."""
    prof = g4hunter_profile(seq)
    if len(prof) <= window:
        return g4hunter_score(seq), 0, len(prof)
    best, at = None, 0
    for i in range(len(prof) - window + 1):
        s = sum(prof[i:i + window]) / window
        if best is None or abs(s) > abs(best):
            best, at = s, i
    return best, at, at + window


def g4_motifs(seq: str):
    seq = seq.upper().replace("U", "T")
    three = sorted({m.group(1) for m in MOTIF_3TETRAD.finditer(seq)})
    two = sorted({m.group(1) for m in MOTIF_2TETRAD.finditer(seq)})
    return {"three_tetrad": three, "two_tetrad": two}


# ---- local alignment (Smith-Waterman) for the artefact control ---------------
def smith_waterman(a: str, b: str, match=2, mismatch=-1, gap=-2):
    """Return (score, aligned_a, aligned_b). Small sequences, so the plain
    O(nm) implementation is fine and avoids a dependency."""
    a, b = a.upper().replace("U", "T"), b.upper().replace("U", "T")
    n, m = len(a), len(b)
    H = [[0] * (m + 1) for _ in range(n + 1)]
    ptr = [[0] * (m + 1) for _ in range(n + 1)]
    best, bi, bj = 0, 0, 0
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = H[i - 1][j - 1] + (match if a[i - 1] == b[j - 1] else mismatch)
            up = H[i - 1][j] + gap
            left = H[i][j - 1] + gap
            val = max(0, diag, up, left)
            H[i][j] = val
            ptr[i][j] = 0 if val == 0 else (1 if val == diag else
                                            (2 if val == up else 3))
            if val > best:
                best, bi, bj = val, i, j
    ra, rb = [], []
    i, j = bi, bj
    while i > 0 and j > 0 and ptr[i][j] != 0:
        d = ptr[i][j]
        if d == 1:
            ra.append(a[i - 1]); rb.append(b[j - 1]); i -= 1; j -= 1
        elif d == 2:
            ra.append(a[i - 1]); rb.append("-"); i -= 1
        else:
            ra.append("-"); rb.append(b[j - 1]); j -= 1
    return best, "".join(reversed(ra)), "".join(reversed(rb))


def alignment_identity(a: str, b: str):
    score, aa, bb = smith_waterman(a, b)
    if not aa:
        return {"score": 0, "aligned_len": 0, "identity": 0.0,
                "query": "", "subject": ""}
    ident = sum(1 for x, y in zip(aa, bb) if x == y and x != "-")
    return {"score": score, "aligned_len": len(aa),
            "identity": round(ident / len(aa), 3),
            "identity_over_ref": round(ident / len(b), 3),
            "query": aa, "subject": bb}


def alignment_pvalue(query: str, ref: str, n: int = 2000, seed: int = 0):
    """How surprising is this alignment score, given the query's composition?

    G-rich sequences align to G-rich sequences by construction, so a high raw
    identity against a G4 aptamer proves nothing on its own. The null here is
    the query's own mononucleotide composition: shuffle the query n times and
    ask how often a shuffle beats the real alignment score. A dinucleotide
    shuffle would be stricter still, but at these lengths it barely moves and
    the mononucleotide null is already the conservative direction for a
    G-rich query (shuffles keep every G, so they align well too).
    """
    import random
    rng = random.Random(seed)
    obs, _, _ = smith_waterman(query, ref)
    chars = list(query.upper().replace("U", "T"))
    ge = 0
    for _ in range(n):
        rng.shuffle(chars)
        s, _, _ = smith_waterman("".join(chars), ref)
        if s >= obs:
            ge += 1
    return {"observed_score": obs, "n_shuffles": n,
            "n_ge_observed": ge, "p_value": round((ge + 1) / (n + 1), 4)}


def analyse(name: str, seq: str, chemistry: str = ""):
    win_score, w0, w1 = g4hunter_windows(seq)
    motifs = g4_motifs(seq)
    refs = {r["name"]: alignment_identity(seq, r["sequence"])
            for r in REFERENCE_G4_APTAMERS}
    forms_g4 = (abs(win_score) >= G4HUNTER_THRESHOLD
                or bool(motifs["three_tetrad"] or motifs["two_tetrad"]))
    return {
        "name": name, "chemistry": chemistry, "length": len(seq),
        "sequence": seq,
        "g4hunter_whole": round(g4hunter_score(seq), 3),
        "g4hunter_windowed_max": round(win_score, 3),
        "g4hunter_window": [w0, w1],
        "motif_3tetrad": motifs["three_tetrad"],
        "motif_2tetrad": motifs["two_tetrad"],
        "predicted_g4": forms_g4,
        "vs_reference_g4_aptamers": refs,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seq", action="append", default=[],
                    help="extra sequence to analyse (repeatable); "
                         "use NAME=SEQUENCE to label it")
    ap.add_argument("--out", default="results/g4_analysis.json")
    args = ap.parse_args()

    rows = []
    for a in KNOWN_APTAMERS:
        if a["sequence"]:
            rows.append(analyse(a["name"], a["sequence"], a["chemistry"]))
    for r in REFERENCE_G4_APTAMERS:
        rows.append(analyse(f"[ref] {r['name']}", r["sequence"],
                            f"reference, target={r['target']}"))
    for s in args.seq:
        name, _, sq = s.partition("=")
        rows.append(analyse(name if sq else "user", sq or name))

    print(f"{'aptamer':<24} {'len':>4} {'G4H':>7} {'G4H_win':>8} "
          f"{'3-tet':>6} {'2-tet':>6}  G4?")
    print("-" * 74)
    for r in rows:
        print(f"{r['name']:<24} {r['length']:>4} {r['g4hunter_whole']:>7.2f} "
              f"{r['g4hunter_windowed_max']:>8.2f} "
              f"{len(r['motif_3tetrad']):>6} {len(r['motif_2tetrad']):>6}  "
              f"{'YES' if r['predicted_g4'] else 'no'}")

    print("\n=== similarity to canonical G4 aptamers (artefact control) ===")
    for r in rows:
        if r["name"].startswith("[ref]"):
            continue
        for ref, al in r["vs_reference_g4_aptamers"].items():
            if al["aligned_len"] >= 8:
                refseq = next(x["sequence"] for x in REFERENCE_G4_APTAMERS
                              if x["name"] == ref)
                pv = alignment_pvalue(r["sequence"], refseq)
                al["null"] = pv
                print(f"\n{r['name']}  vs  {ref}")
                print(f"  aligned {al['aligned_len']} nt, "
                      f"identity {al['identity']:.0%} "
                      f"({al['identity_over_ref']:.0%} of the reference)")
                print(f"    query   {al['query']}")
                bars = "".join("|" if x == y and x != "-" else " "
                               for x, y in zip(al["query"], al["subject"]))
                print(f"            {bars}")
                print(f"    ref     {al['subject']}")
                print(f"  composition-shuffled null: p = {pv['p_value']} "
                      f"({pv['n_ge_observed']}/{pv['n_shuffles']} shuffles "
                      f"reach score {pv['observed_score']})")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rows, indent=2))
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
