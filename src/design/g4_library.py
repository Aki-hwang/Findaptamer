"""
Focused MMP-9 aptamer library: a computable search-space reduction.

The honest constraint (docs/06, docs/07): we cannot compute "does this
sequence bind MMP-9". The co-folding oracle scored a known 20 nM binder below
random sequence, its pose had zero seed-to-seed overlap, and the same pipeline
scored a canonical protein-DNA complex at ipTM 0.974 +/- 0.0007 -- so the
failure is real and specific to aptamers, not ours to fix by tuning.

So this module never predicts binding. It does the one thing that IS
computable and that a wet SELEX round cannot do cheaply: shrink the search
space from 4^30 (~10^18) to ~10^3-10^4 sequences, using only quantities that
are measurable or algorithmically defined.

The reduction rests on one experimental observation, not on a model: three
independent selections against MMP-9 converged on G-quadruplex folds (8F14A
confirmed by CD/Tm, the LVMH patent series claimed as G4, and MMP9-DNA-30
scoring +1.44 by G4Hunter). So the scaffold is fixed to a G4 topology and the
LOOPS are varied -- the loops are the recognition surface (in TBA the TT
loops, not the G-tracts, contact thrombin), and they are exactly what a
selection would explore.

That makes this the computational half of ML-guided SELEX: this library goes
into one wet round, the round data trains a RaptRanker/DeepAptamer-class model,
and the model then runs the in-silico rounds. The literature's working recipe
(DL-SELEX: up to 450x affinity, 80% fewer rounds) needs target-specific SELEX
data, and this is how we generate it with 10^4 sequences instead of 10^15.

Filters, in the order they are applied (cheapest first):
  1. G4 topology       -- 2-tetrad motif must be present (regex, deterministic)
  2. flank penalty     -- G4Hunter (Bedrat & Mergny 2016) is REPORTED but not
                          gated on the core: within a scaffold-fixed library
                          with G-free loops it reduces to 16/length, a
                          loop-length filter in disguise. It is gated only on
                          how much the flanks degrade the core, which is what
                          it can actually measure here. See screen().
  3. competing fold    -- ViennaRNA DNA params: a hairpin too stable to let the
                          G4 win is disqualifying. G4s are invisible to the
                          nearest-neighbour model, so a very negative MFE here
                          means a DUPLEX competitor, not a good aptamer.
  4. TBA novelty       -- composition-preserving permutation test. A candidate
                          that is just the thrombin aptamer must not enter the
                          library: MMP9-DNA-30 already carries that liability
                          (p = 0.0035, docs/07) and its assays run in tears and
                          serum where thrombin is present.
  5. synthesis sanity  -- no run of >4 identical bases outside the G-tracts,
                          which are hard to synthesise and prone to slippage.
  6. diversity         -- greedy max-min Hamming selection over the survivors,
                          so the library spans loop space instead of clustering.

Usage:
    python src/design/g4_library.py --n 2000 --out results/library_v1.csv
    python src/design/g4_library.py --n 200 --tetrads 3
"""
from __future__ import annotations

import argparse
import csv
import json
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from analysis.g4 import (  # noqa: E402
    alignment_identity, alignment_pvalue, g4_motifs, g4hunter_windows,
)
from target.mmp9_aptamers import REFERENCE_G4_APTAMERS  # noqa: E402

TBA = next(r["sequence"] for r in REFERENCE_G4_APTAMERS
           if r["name"].startswith("TBA"))

# Loop-length ranges for a chair-type intramolecular G4. Lateral loops are
# short (2-3 nt in TBA); the central loop spans the groove and tolerates more.
LATERAL_LOOP_LEN = (2, 3)
CENTRAL_LOOP_LEN = (3, 5)
# Loops are drawn without G to keep the G-tract register unambiguous: a G in a
# loop can slip into a tetrad and change the fold, which would silently break
# the one structural assumption the whole library rests on.
LOOP_ALPHABET = "ACT"
FLANK_ALPHABET = "ACGT"


def _rand(rng, alphabet, n):
    return "".join(rng.choice(alphabet) for _ in range(n))


def make_candidate(rng, tetrads: int, flank5: int, flank3: int) -> dict:
    """Build one G4 scaffold with randomised loops and flanks."""
    tract = "G" * tetrads
    l1 = _rand(rng, LOOP_ALPHABET, rng.randint(*LATERAL_LOOP_LEN))
    l2 = _rand(rng, LOOP_ALPHABET, rng.randint(*CENTRAL_LOOP_LEN))
    l3 = _rand(rng, LOOP_ALPHABET, rng.randint(*LATERAL_LOOP_LEN))
    f5 = _rand(rng, FLANK_ALPHABET, flank5) if flank5 else ""
    f3 = _rand(rng, FLANK_ALPHABET, flank3) if flank3 else ""
    core = f"{tract}{l1}{tract}{l2}{tract}{l3}{tract}"
    return {"sequence": f5 + core + f3, "core": core,
            "loop1": l1, "loop2": l2, "loop3": l3,
            "flank5": f5, "flank3": f3, "tetrads": tetrads}


# ---- filters -----------------------------------------------------------------
def competing_fold_mfe(seq: str):
    """MFE under DNA nearest-neighbour parameters, or None if unavailable.

    Interpretation is inverted relative to normal RNA design. ViennaRNA has no
    G-quadruplex term, so it can only see Watson-Crick competitors. A near-zero
    MFE therefore means "nothing competes with the G4" -- which is what we want.
    A strongly negative MFE means a stem-loop that would sequester the G-tracts.
    """
    try:
        import RNA
    except ImportError:
        return None
    RNA.params_load_DNA_Mathews2004()
    fc = RNA.fold_compound(seq, RNA.md())
    _, mfe = fc.mfe()
    return round(mfe, 2)


def homopolymer_ok(cand: dict, max_run: int = 4) -> bool:
    """No long identical runs outside the deliberate G-tracts."""
    seq, i = cand["sequence"], 0
    while i < len(seq):
        j = i
        while j < len(seq) and seq[j] == seq[i]:
            j += 1
        run = j - i
        if run > max_run and not (seq[i] == "G" and run <= cand["tetrads"]):
            return False
        i = j
    return True


def screen(cand: dict, g4_min: float, mfe_floor: float, tba_p_min: float,
           n_shuffles: int) -> tuple[bool, str, dict]:
    """Apply the filters cheapest-first.

    Returns (passed, code, metrics) where `code` is a FIXED string so the
    caller can tally rejections by filter. Embedding the measured value in the
    code would give every rejection its own bucket and hide the distribution.
    """
    seq = cand["sequence"]
    m = {}

    motifs = g4_motifs(seq)
    m["motif_2tetrad"] = len(motifs["two_tetrad"])
    m["motif_3tetrad"] = len(motifs["three_tetrad"])
    if not (motifs["two_tetrad"] or motifs["three_tetrad"]):
        return False, "no G4 motif", m

    # G4Hunter is REPORTED, not gated, for the core -- and the reason matters.
    # Every candidate here is built on the same scaffold with G-free loops, so
    # its G4Hunter score is exactly (4 * tetrads * 2) / length: a monotone
    # function of loop length carrying no information about G4 propensity
    # within this library. TBA scores 1.133 only because it happens to have a
    # lone G in its TGT central loop. Gating on it would quietly select for
    # short loops -- a bias dressed up as a structural criterion. The topology
    # regex above already enforces what G4Hunter was standing in for.
    core_win, _, _ = g4hunter_windows(cand["core"])
    full_win, _, _ = g4hunter_windows(seq)
    m["g4hunter_core"] = round(core_win, 3)
    m["g4hunter"] = round(full_win, 3)
    m["flank_penalty"] = round(core_win - full_win, 3)

    # The one part of G4Hunter that IS length-fair: its SIGN. A negative core
    # means the loops carry more C character than the tracts carry G character,
    # and loop C's are the complement of the tracts -- they can pair with them
    # and dismantle the fold. Unlike the magnitude, the sign does not move with
    # length, so gating on it smuggles in no loop-length preference.
    if core_win <= 0:
        return False, "C-rich loops outweigh the G-tracts", m
    # Belt and braces, mechanistically rather than statistically: a CCC+ run
    # inside a loop is a direct pairing liability even when the overall sign
    # survives, and it is also an i-motif seed at low pH.
    for loop in (cand["loop1"], cand["loop2"], cand["loop3"]):
        if "CCC" in loop:
            return False, "CCC run inside a loop", m
    # Where G4Hunter IS informative: flanks are drawn from the full alphabet, so
    # a C-rich flank can both score negative and pair with a G-tract. Gate on
    # the degradation the flanks cause, which is length-fair by construction.
    if (cand["flank5"] or cand["flank3"]) and core_win - full_win > g4_min:
        return False, "flanks degrade the G4 core", m

    mfe = competing_fold_mfe(seq)
    m["competing_mfe"] = mfe
    if mfe is not None and mfe < mfe_floor:
        return False, "competing Watson-Crick hairpin", m

    if not homopolymer_ok(cand):
        return False, "homopolymer run > 4", m

    # Most expensive filter last: only sequences that already look like usable
    # G4s pay for hundreds of alignments.
    al = alignment_identity(seq, TBA)
    m["tba_identity"] = al["identity"]
    m["tba_aligned_len"] = al["aligned_len"]
    pv = alignment_pvalue(seq, TBA, n=n_shuffles)
    m["tba_p"] = pv["p_value"]
    if pv["p_value"] < tba_p_min:
        return False, "too close to TBA", m

    return True, "pass", m


def hamming_pad(a: str, b: str) -> int:
    """Hamming distance with the shorter sequence padded — candidates differ in
    length, and a length difference is itself a real difference."""
    n = max(len(a), len(b))
    a, b = a.ljust(n, "-"), b.ljust(n, "-")
    return sum(1 for x, y in zip(a, b) if x != y)


def diversify(rows: list[dict], k: int) -> list[dict]:
    """Greedy max-min selection: repeatedly take the candidate farthest from
    everything already chosen. Random sampling would over-represent whatever
    loop compositions are most probable; this spans the space instead."""
    if len(rows) <= k:
        return rows
    chosen = [max(rows, key=lambda r: r["metrics"]["g4hunter"])]
    remaining = [r for r in rows if r is not chosen[0]]
    mind = {id(r): hamming_pad(r["sequence"], chosen[0]["sequence"])
            for r in remaining}
    while len(chosen) < k and remaining:
        best = max(remaining, key=lambda r: mind[id(r)])
        chosen.append(best)
        remaining.remove(best)
        for r in remaining:
            d = hamming_pad(r["sequence"], best["sequence"])
            if d < mind[id(r)]:
                mind[id(r)] = d
    return chosen


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n", type=int, default=1000, help="library size to emit")
    ap.add_argument("--tetrads", type=int, default=2, choices=(2, 3),
                    help="2 = TBA-class chair (the MMP-9 precedent); "
                         "3 = three-tetrad, more stable but larger footprint")
    ap.add_argument("--flank5", type=int, default=0)
    ap.add_argument("--flank3", type=int, default=0)
    ap.add_argument("--oversample", type=int, default=12,
                    help="candidates generated per library slot before filtering")
    ap.add_argument("--max-flank-penalty", type=float, default=0.2,
                    dest="g4_min",
                    help="max G4Hunter drop the flanks may cause "
                         "(core score minus full-sequence score). Only applied "
                         "when flanks are requested; see screen() for why the "
                         "raw G4Hunter score is reported rather than gated")
    ap.add_argument("--mfe-floor", type=float, default=-3.0,
                    help="reject if a Watson-Crick competitor is more stable "
                         "than this (kcal/mol)")
    ap.add_argument("--tba-p-min", type=float, default=0.05,
                    help="reject candidates whose TBA similarity is significant")
    ap.add_argument("--shuffles", type=int, default=500)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="results/g4_library.csv")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    target_pool = args.n * args.oversample
    seen, passed, rejects = set(), [], {}

    print(f"generating {target_pool} {args.tetrads}-tetrad candidates "
          f"-> filtering -> {args.n} diverse survivors", flush=True)
    for i in range(target_pool):
        cand = make_candidate(rng, args.tetrads, args.flank5, args.flank3)
        if cand["sequence"] in seen:
            continue
        seen.add(cand["sequence"])
        ok, reason, metrics = screen(cand, args.g4_min, args.mfe_floor,
                                     args.tba_p_min, args.shuffles)
        if ok:
            passed.append({**cand, "metrics": metrics})
        else:
            rejects[reason] = rejects.get(reason, 0) + 1
        if (i + 1) % 2000 == 0:
            print(f"  {i+1}/{target_pool} generated, {len(passed)} passed",
                  flush=True)

    print(f"\n{len(seen)} unique candidates -> {len(passed)} passed all filters "
          f"({100*len(passed)/max(len(seen),1):.1f}%)")
    print("rejections by first failing filter:")
    for k, v in sorted(rejects.items(), key=lambda x: -x[1]):
        print(f"  {v:>7}  {k}")

    library = diversify(passed, args.n)
    if library:
        d = [hamming_pad(a["sequence"], b["sequence"])
             for i, a in enumerate(library) for b in library[i + 1:i + 6]]
        print(f"\nlibrary: {len(library)} sequences, "
              f"lengths {min(len(r['sequence']) for r in library)}-"
              f"{max(len(r['sequence']) for r in library)} nt, "
              f"mean pairwise Hamming {sum(d)/len(d):.1f}" if d else "")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["id", "sequence", "length", "loop1", "loop2", "loop3",
                    "g4hunter_core", "g4hunter_full", "flank_penalty",
                    "competing_mfe", "tba_identity", "tba_p"])
        for i, r in enumerate(library):
            m = r["metrics"]
            w.writerow([f"MMP9G4_{i:05d}", r["sequence"], len(r["sequence"]),
                        r["loop1"], r["loop2"], r["loop3"],
                        m["g4hunter_core"], m["g4hunter"], m["flank_penalty"],
                        m["competing_mfe"], m["tba_identity"], m["tba_p"]])
    print(f"\nWrote {out}")

    meta = out.with_suffix(".meta.json")
    meta.write_text(json.dumps({
        "n_generated": len(seen), "n_passed": len(passed),
        "n_emitted": len(library), "rejections": rejects,
        "params": vars(args),
        "provenance": "G4 prior from three independent MMP-9 selections; "
                      "see docs/07_g4_discovery_ko.md",
        "what_this_is_not": "This library is NOT a binding prediction. No "
                            "filter here estimates affinity for MMP-9. It is a "
                            "search-space reduction for one wet SELEX round.",
    }, indent=2))
    print(f"Wrote {meta}")


if __name__ == "__main__":
    main()
