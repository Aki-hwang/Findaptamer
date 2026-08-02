"""
Analyse REAL MMP-9 SELEX output — and use it to falsify our own assumptions.

Everything downstream in this project rests on one inherited hypothesis: that
MMP-9 aptamers are G-quadruplexes. That came from three published selections
(docs/07), but two of them share a laboratory, and the third piece of evidence
is our own G4Hunter prediction rather than an experiment. It has never been
tested against an independent selection.

Experimental SELEX output against MMP-9 tests it directly, so this script is
written to let the hypothesis LOSE. It reports the G4 fraction of the real hits
against a composition-matched null built from the same sequences, and states
plainly when the prior is not supported. A G4 prior that survives this is worth
something; one that is only ever confirmed by the analysis that assumes it is
not.

Four things it checks:

  1. G4 prior          what fraction of real hits form G4, versus a
                       dinucleotide-shuffled null of the same sequences. Shuffle
                       preserves base and neighbour composition, so a G-rich
                       library does not manufacture a positive on its own.
  2. TBA artefact      the thrombin-aptamer permutation test from docs/07.
                       MMP9-DNA-30, the published DNA binder, contains TBA at
                       p=0.0035. If real hits do too, that is a selection
                       artefact to catch BEFORE anyone characterises them.
  3. known binders     alignment to MMP9-DNA-30 and F3B -- did this selection
                       rediscover the published aptamers or find new chemistry?
  4. our predictions   whether the focused library (results/g4_library_*.csv)
                       and the blind patent prediction overlap the real hits.
                       This scores our design retrospectively; a miss is
                       reportable and is not quietly dropped.

Input: FASTA, or CSV/TSV with a sequence column and optionally a count and a
round column. Primer constant regions are stripped when --fwd/--rev are given,
because a constant region is in every read and would dominate every statistic.

Usage:
    python src/analysis/selex_hits.py --seqs hits.csv
    python src/analysis/selex_hits.py --seqs hits.fasta \
        --fwd ATTCCAACCTGCAAGAACAA --rev TTACTAAGCGAGACGTCGTA
"""
from __future__ import annotations

import argparse
import csv
import json
import random
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from analysis.g4 import (  # noqa: E402
    alignment_identity, alignment_pvalue, g4_motifs, g4hunter_windows,
)
from target.mmp9_aptamers import KNOWN_APTAMERS, REFERENCE_G4_APTAMERS  # noqa: E402

COMP = str.maketrans("ACGT", "TGCA")
TBA = next(r["sequence"] for r in REFERENCE_G4_APTAMERS
           if r["name"].startswith("TBA"))


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def load_sequences(path: str):
    """FASTA or delimited text. Returns [{sequence, count, round, name}]."""
    p = Path(path)
    text = p.read_text()
    rows = []
    if text.lstrip().startswith(">"):
        name, buf = None, []
        for line in text.splitlines():
            if line.startswith(">"):
                if name is not None:
                    rows.append({"name": name, "sequence": "".join(buf)})
                name, buf = line[1:].strip(), []
            elif line.strip():
                buf.append(line.strip())
        if name is not None:
            rows.append({"name": name, "sequence": "".join(buf)})
        for i, r in enumerate(rows):
            r.setdefault("count", None)
            r.setdefault("round", None)
            # a "seq_count" style header is common from NGS clustering
            m = re.search(r"(?:count|reads?|freq)[=_:\s]+(\d+)", r["name"], re.I)
            if m:
                r["count"] = int(m.group(1))
        return rows

    delim = "\t" if "\t" in text.splitlines()[0] else ","
    rd = csv.DictReader(text.splitlines(), delimiter=delim)
    def pick(d, *names):
        for k in d:
            if k and k.strip().lower() in names:
                return d[k]
        return None
    for i, d in enumerate(rd):
        s = pick(d, "sequence", "seq", "aptamer", "dna")
        if not s:
            continue
        c = pick(d, "count", "reads", "frequency", "freq", "n")
        rnd = pick(d, "round", "cycle", "r")
        rows.append({"name": pick(d, "name", "id") or f"seq{i}",
                     "sequence": s.strip().upper(),
                     "count": int(c) if c and str(c).strip().isdigit() else None,
                     "round": rnd})
    return rows


def strip_constant(seq: str, fwd: str | None, rev: str | None):
    """Remove the 5' and 3' constant regions if present (either orientation)."""
    s = seq.upper().replace("U", "T")
    if fwd and s.startswith(fwd.upper()):
        s = s[len(fwd):]
    if rev:
        tail = revcomp(rev.upper())
        if s.endswith(tail):
            s = s[:-len(tail)]
        elif s.endswith(rev.upper()):
            s = s[:-len(rev)]
    return s


def dinuc_shuffle(seq: str, rng: random.Random) -> str:
    """Altschul-Erikson dinucleotide shuffle.

    A plain mononucleotide shuffle is too weak a null here: G-quadruplex
    propensity depends on G RUNS, and preserving only base composition breaks
    them up, so almost any G-rich sequence would look significantly G4 against
    it. Preserving dinucleotide counts keeps run structure and makes the test
    honest.
    """
    if len(seq) < 4:
        return seq
    edges = {}
    for a, b in zip(seq, seq[1:]):
        edges.setdefault(a, []).append(b)
    for k in edges:
        rng.shuffle(edges[k])
    # walk an Eulerian path; fall back to the original on a dead end
    out, cur = [seq[0]], seq[0]
    pool = {k: list(v) for k, v in edges.items()}
    for _ in range(len(seq) - 1):
        if not pool.get(cur):
            return seq
        nxt = pool[cur].pop()
        out.append(nxt)
        cur = nxt
    return "".join(out)


def is_g4(seq: str) -> bool:
    m = g4_motifs(seq)
    score, _, _ = g4hunter_windows(seq)
    return bool(m["two_tetrad"] or m["three_tetrad"]) and score > 0


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seqs", required=True)
    ap.add_argument("--fwd", help="5' constant region to strip")
    ap.add_argument("--rev", help="reverse primer (its revcomp is stripped)")
    ap.add_argument("--top", type=int, default=25)
    ap.add_argument("--shuffles", type=int, default=500)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--library", default="results/g4_library_v1.csv")
    ap.add_argument("--out", default="results/selex_hits_analysis.json")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    rows = load_sequences(args.seqs)
    if not rows:
        raise SystemExit(f"no sequences parsed from {args.seqs}")
    for r in rows:
        r["core"] = strip_constant(r["sequence"], args.fwd, args.rev)
    rows.sort(key=lambda r: -(r["count"] or 0))

    n_stripped = sum(1 for r in rows if len(r["core"]) < len(r["sequence"]))
    print(f"{len(rows)} sequences | constant region stripped from {n_stripped}")
    lens = [len(r["core"]) for r in rows]
    print(f"core length: min {min(lens)} med {sorted(lens)[len(lens)//2]} "
          f"max {max(lens)}\n")

    # ---- 1. does the G4 prior survive? ------------------------------------
    obs = sum(1 for r in rows if is_g4(r["core"]))
    null = []
    for _ in range(max(args.shuffles // 50, 20)):
        null.append(sum(1 for r in rows if is_g4(dinuc_shuffle(r["core"], rng))))
    null_mean = sum(null) / len(null)
    ge = sum(1 for v in null if v >= obs)
    p_g4 = (ge + 1) / (len(null) + 1)
    frac = obs / len(rows)
    print("=== 1. G4 prior ===")
    print(f"  real hits forming G4      : {obs}/{len(rows)}  ({frac:.0%})")
    print(f"  dinucleotide-shuffled null: {null_mean:.1f}/{len(rows)}  "
          f"({null_mean/len(rows):.0%})")
    print(f"  p = {p_g4:.4f}")
    g4_supported = p_g4 < 0.05 and frac > null_mean / len(rows)
    print(f"  -> G4 prior {'SUPPORTED' if g4_supported else 'NOT SUPPORTED'}"
          f" by this selection")
    if not g4_supported:
        print("     The prior came from published selections, two of which share")
        print("     a lab. This independent selection does not back it, so the")
        print("     G4-focused library design does not apply to these hits.")

    # ---- 2. TBA artefact ---------------------------------------------------
    print("\n=== 2. thrombin-aptamer artefact check ===")
    tba_hits = []
    for r in rows[:args.top]:
        al = alignment_identity(r["core"], TBA)
        if al["aligned_len"] >= 8:
            pv = alignment_pvalue(r["core"], TBA, n=args.shuffles)
            r["tba_p"] = pv["p_value"]
            r["tba_identity"] = al["identity"]
            if pv["p_value"] < 0.05:
                tba_hits.append(r)
    if tba_hits:
        print(f"  {len(tba_hits)} of the top {min(args.top,len(rows))} match TBA "
              f"significantly:")
        for r in tba_hits[:5]:
            print(f"    {r['name']:<20} identity {r['tba_identity']:.0%}  "
                  f"p={r['tba_p']}")
        print("  -> counter-screen these against thrombin before characterising.")
    else:
        print(f"  none of the top {min(args.top,len(rows))} match TBA "
              f"significantly (good)")

    # ---- 3. published MMP-9 binders ---------------------------------------
    print("\n=== 3. vs published MMP-9 aptamers ===")
    known = [(a["name"], a["sequence"]) for a in KNOWN_APTAMERS if a["sequence"]]
    best_known = []
    for name, kseq in known:
        best = max(rows[:args.top],
                   key=lambda r: alignment_identity(r["core"], kseq)["score"])
        al = alignment_identity(best["core"], kseq)
        best_known.append({"known": name, "closest_hit": best["name"],
                           "identity": al["identity"],
                           "aligned_len": al["aligned_len"]})
        print(f"  {name:<14} closest hit {best['name']:<18} "
              f"identity {al['identity']:.0%} over {al['aligned_len']} nt")

    # ---- 4. did OUR design anticipate these? ------------------------------
    print("\n=== 4. our designed library vs the real hits ===")
    lib_path = Path(args.library)
    overlap = None
    if lib_path.exists():
        lib = [r["sequence"] for r in csv.DictReader(lib_path.read_text()
                                                     .splitlines())]
        exact = {r["core"] for r in rows} & set(lib)
        best_sim = 0.0
        for r in rows[:args.top]:
            for L in lib[:400]:
                al = alignment_identity(r["core"], L)
                if al["aligned_len"] >= 10:
                    best_sim = max(best_sim, al["identity"])
        overlap = {"exact_matches": len(exact),
                   "best_identity_to_library": round(best_sim, 3)}
        print(f"  exact matches to results/g4_library_v1.csv: {len(exact)}")
        print(f"  best identity of any real hit to the library: {best_sim:.0%}")
        if len(exact) == 0:
            print("  -> our library did not contain these sequences. Expected "
                  "for a 10^24 space, but it means the design is unvalidated.")
    else:
        print(f"  {lib_path} not found — skipped")

    payload = {
        "n_sequences": len(rows),
        "n_constant_stripped": n_stripped,
        "g4_prior": {"observed_g4": obs, "n": len(rows),
                     "fraction": round(frac, 3),
                     "null_mean": round(null_mean, 2),
                     "p_value": round(p_g4, 4),
                     "supported": bool(g4_supported)},
        "tba_artefact_hits": [{"name": r["name"], "p": r.get("tba_p"),
                               "identity": r.get("tba_identity")}
                              for r in tba_hits],
        "vs_published": best_known,
        "vs_our_library": overlap,
        "top_hits": [{"name": r["name"], "count": r["count"],
                      "core": r["core"], "length": len(r["core"]),
                      "g4hunter": round(g4hunter_windows(r["core"])[0], 3),
                      "is_g4": is_g4(r["core"])}
                     for r in rows[:args.top]],
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2))
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
