"""
Fragment pool construction for the assembler.

AiDTA builds its 44-fragment pool by docking all 10,880 fragments to the target
epitope with HDOCK and keeping the top 5 ss + top 5 ds per length (+ the 4
mononucleotides). That docking step exists because AiDTA's reward never sees the
protein — the pool is the *only* place target information enters.

Our pipeline puts the protein in the reward (Boltz-2 co-folding oracle), so we
support two modes:

  build_generic_pool()  — target-agnostic, diversity-maximizing pool. Same shape
                          as AiDTA's (44 fragments) but selected for sequence
                          diversity instead of docking score. Runs instantly, no
                          docking needed; the oracle does the targeting.

  build_docked_pool()   — AiDTA-style: read HDOCK results and keep the best
                          epitope-binding fragments. Use this when docking output
                          is available, to seed the search with target bias.

Both return a list of (fragment_sequence, fragment_structure) usable directly by
src/generator/assembler.py.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from fragments.build_library import enumerate_ss, enumerate_ds  # noqa: E402

MONONUCLEOTIDES = [("A", "."), ("G", "."), ("C", "."), ("T", ".")]


def _hamming(a: str, b: str) -> int:
    return sum(1 for x, y in zip(a, b) if x != y)


def _max_min_diverse(items, k: int):
    """Greedy max-min selection: repeatedly take the item farthest (Hamming) from
    everything already chosen. Deterministic — seeds with the first item.

    `items` is a list of (sequence, structure); comparison uses the first strand
    so that ds fragments are compared on their defining strand.
    """
    if k >= len(items):
        return list(items)

    def key(seq):
        return seq.split("&")[0]

    chosen = [items[0]]
    remaining = list(items[1:])
    while len(chosen) < k and remaining:
        best, best_d = None, -1
        for cand in remaining:
            d = min(_hamming(key(cand[0]), key(c[0])) for c in chosen)
            if d > best_d:
                best, best_d = cand, d
        chosen.append(best)
        remaining.remove(best)
    return chosen


def build_generic_pool(n_per_class: int = 5, lengths=(3, 4, 5, 6),
                       include_mono: bool = True):
    """Target-agnostic pool with AiDTA's shape (default: 5 ss + 5 ds per length
    + 4 mononucleotides = 44 fragments), selected to maximize sequence diversity.
    """
    pool = []
    for L in lengths:
        ss = [f for f in enumerate_ss((L,))]
        ds = [f for f in enumerate_ds((L,))]
        pool.extend(_max_min_diverse(ss, n_per_class))
        pool.extend(_max_min_diverse(ds, n_per_class))
    if include_mono:
        pool = MONONUCLEOTIDES + pool
    return pool


def build_docked_pool(docking_csv: str, n_per_class: int = 5,
                      lengths=(3, 4, 5, 6), include_mono: bool = True,
                      score_col: str = "docking_score",
                      ascending: bool = True):
    """AiDTA-style pool from docking results.

    `docking_csv` must have columns: sequence, structure, type ('ss'|'ds'),
    length (nt or bp), and a score column (default 'docking_score'). HDOCK scores
    are more negative = better, so ascending=True keeps the best.
    """
    rows = []
    with open(docking_csv) as f:
        for r in csv.DictReader(f):
            try:
                r["_score"] = float(r[score_col])
                r["_len"] = int(r.get("length") or r.get("bp"))
            except (TypeError, ValueError, KeyError):
                continue
            rows.append(r)

    pool = []
    for L in lengths:
        for typ in ("ss", "ds"):
            sub = [r for r in rows if r.get("type") == typ and r["_len"] == L]
            sub.sort(key=lambda r: r["_score"], reverse=not ascending)
            for r in sub[:n_per_class]:
                pool.append((r["sequence"], r["structure"]))
    if include_mono:
        pool = MONONUCLEOTIDES + pool
    return pool


def save_pool(pool, path: str):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sequence", "structure"])
        w.writerows(pool)


def load_pool(path: str):
    with open(path) as f:
        return [(r["sequence"], r["structure"]) for r in csv.DictReader(f)]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mode", choices=["generic", "docked"], default="generic")
    ap.add_argument("--docking-csv", help="required for --mode docked")
    ap.add_argument("--n-per-class", type=int, default=5)
    ap.add_argument("--out", default="data/mmp9/pool.csv")
    args = ap.parse_args()

    if args.mode == "docked":
        if not args.docking_csv:
            ap.error("--docking-csv is required for --mode docked")
        pool = build_docked_pool(args.docking_csv, args.n_per_class)
    else:
        pool = build_generic_pool(args.n_per_class)

    save_pool(pool, args.out)
    ss = sum(1 for s, _ in pool if "&" not in s and len(s) > 1)
    ds = sum(1 for s, _ in pool if "&" in s)
    mono = sum(1 for s, _ in pool if len(s) == 1)
    print(json.dumps({"mode": args.mode, "total": len(pool),
                      "mononucleotides": mono, "ss": ss, "ds": ds,
                      "out": args.out}, indent=2))


if __name__ == "__main__":
    main()
