"""
Three MMP-9 aptamer candidates, pre-registered before seeing the user's SELEX.

This file is the timestamp. The candidates, the reasoning behind each, and the
grading rule are all fixed here and committed BEFORE the experimental sequences
are shown to us. Scoring a prediction against a rule invented after seeing the
answer is worthless, and this project has already used that discipline once for
the LVMH patent (src/predict/lvmh_patent.py).

WHAT WE CANNOT DO, stated up front. There is no validated binding oracle for
MMP-9. Four independent measurements, each in this repository:

  co-folding ipTM        two known binders both scored BELOW random negatives;
                         signal/noise ~1/40 while the same pipeline scored a
                         canonical protein-DNA complex at 0.974 +/- 0.0007
                         (docs/06, jobs 3406 / 3410 / 3408)
  co-folding pose        positive within-seed contact Jaccard 0.0
  sequence ML            protein-grouped AUROC 0.210 -- below chance, because
                         each aptamer is positive for exactly one protein, so a
                         novel target has no association to recall (docs/08)
  Apta-MCTS scorers      AUROC 0.4927 and 0.5055 on real pairs, i.e. chance,
                         measured directly on the pretrained weights

So none of the three candidates below carries a predicted affinity, because any
number we printed would be fabricated. They come from a STRUCTURAL PRIOR and are
labelled with the arm that produced them.

THE PRIOR, and its weakness. Every published MMP-9 aptamer whose fold is known
is a G-quadruplex: 8F14A (G-quartet by CD/Tm), the LVMH inhibitor patent (G4
claimed), and MMP9-DNA-30 (G4Hunter +1.44). The counter-example F3B is not, and
scores -1.17, which is what makes this a real signal rather than a tautology.
The weakness is that 8F14A and the LVMH series come from the SAME laboratory
(Toulme/Dausse), so shared library or protocol bias is a live alternative
explanation (docs/07).

THE HEDGE. Because the prior could be wrong in two different ways, the three
candidates deliberately span three hypotheses rather than being three variants
of one guess:

  A  two-tetrad chair G4   the topology of MMP9-DNA-30's validated core, which
                           is the only unmodified-DNA MMP-9 binder reused by
                           three independent groups
  B  parallel (GGGT)n G4   the other major G4 topology, the class the LVMH
                           inhibitor patent claims, and the class our
                           leave-one-out-validated predictor favours
  C  machine learning      Apta-MCTS output against full-length MMP-9, kept
                           even though its scorer measured at chance, because
                           the request was explicitly for an ML arm and an
                           honest comparison needs it

Both A and B must pass the thrombin-aptamer novelty test (docs/07): MMP9-DNA-30
contains TBA at p = 0.0035, and a candidate that is merely the thrombin aptamer
would cross-react in serum and tears.

Usage:
    python src/predict/mmp9_top3.py --mcts <path>/MMP9_*.csv
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
from target.mmp9_aptamers import KNOWN_APTAMERS, REFERENCE_G4_APTAMERS  # noqa: E402

TBA = next(r["sequence"] for r in REFERENCE_G4_APTAMERS
           if r["name"].startswith("TBA"))
DNA30 = next(a["sequence"] for a in KNOWN_APTAMERS if a["name"] == "MMP9-DNA-30")
TARGET_LEN = 40           # match an N40 selection
LOOP_ALPHABET = "ACT"     # G excluded so the tract register stays unambiguous


def competing_mfe(seq: str):
    try:
        import RNA
    except ImportError:
        return None
    RNA.params_load_DNA_Mathews2004()
    return round(RNA.fold_compound(seq, RNA.md()).mfe()[1], 2)


def build_core(rng, tetrads: int, loops: tuple[int, int, int]) -> dict:
    t = "G" * tetrads
    l1 = "".join(rng.choice(LOOP_ALPHABET) for _ in range(loops[0]))
    l2 = "".join(rng.choice(LOOP_ALPHABET) for _ in range(loops[1]))
    l3 = "".join(rng.choice(LOOP_ALPHABET) for _ in range(loops[2]))
    return {"core": f"{t}{l1}{t}{l2}{t}{l3}{t}", "loops": (l1, l2, l3)}


def pad_to(seq: str, n: int, rng) -> str:
    """Extend to n nt with short G-free flanks.

    Flanks are kept SHORT on purpose. The first version of this design padded a
    15-nt G4 core out to 40 nt with 25 nt of A/T filler, which is not an aptamer
    -- it is a small motif floating in junk. In a real selection the loops are
    the recognition surface (in TBA it is the TT loops, not the G-tracts, that
    contact thrombin), so a 40-mer should spend its length on loops rather than
    on padding. build_core now does that and this only tops up the remainder.

    Flank composition: no G, so a flank cannot donate a tetrad column and shift
    the register; A/T/C only, and C is capped by the caller's competing-fold
    check rather than banned, so the flanks are not a monotonous A/T tract.
    """
    need = n - len(seq)
    if need <= 0:
        return seq[:n]
    left = need // 2
    alpha = "ACT"
    f5 = "".join(rng.choice(alpha) for _ in range(left))
    f3 = "".join(rng.choice(alpha) for _ in range(need - left))
    return f5 + seq + f3


def score_candidate(seq: str, shuffles: int = 1000) -> dict:
    g4h, _, _ = g4hunter_windows(seq)
    m = g4_motifs(seq)
    tba = alignment_identity(seq, TBA)
    tba_p = alignment_pvalue(seq, TBA, n=shuffles)["p_value"]
    d30 = alignment_identity(seq, DNA30)
    return {
        "sequence": seq, "length": len(seq),
        "g4hunter": round(g4h, 3),
        "motif_2tetrad": len(m["two_tetrad"]),
        "motif_3tetrad": len(m["three_tetrad"]),
        "competing_mfe": competing_mfe(seq),
        "tba_identity": tba["identity"], "tba_p": tba_p,
        "vs_MMP9_DNA_30_identity": d30["identity"],
    }


def generate_arm(rng, tetrads, loops_choices, n_try=3000, shuffles=400):
    """Best candidate of a topology: G4-positive, novel vs TBA, no competitor."""
    best = None
    for _ in range(n_try):
        loops = tuple(rng.choice(c) for c in loops_choices)
        c = build_core(rng, tetrads, loops)
        seq = pad_to(c["core"], TARGET_LEN, rng)
        if not (g4_motifs(seq)["two_tetrad"] or g4_motifs(seq)["three_tetrad"]):
            continue
        mfe = competing_mfe(seq)
        if mfe is not None and mfe < -2.0:      # a stable duplex would win
            continue
        pv = alignment_pvalue(seq, TBA, n=shuffles)["p_value"]
        if pv < 0.05:                            # thrombin-aptamer artefact
            continue
        g4h, _, _ = g4hunter_windows(seq)
        if best is None or g4h > best[0]:
            best = (g4h, seq, c["loops"])
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mcts", help="Apta-MCTS output CSV (the ML arm)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="results/mmp9_top3_prediction.json")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    cands = []

    # Loop lengths are chosen so the G4 core FILLS most of the 40 nt. Two
    # tetrads give 8 G, so loops of 6-7 nt each put the core at 26-29 nt and
    # leave only ~11 nt of flank instead of the 25 nt of filler the first
    # version produced. Long loops are also the point biologically: they are
    # the surface that contacts protein, and a 40-mer can afford them.
    print("=== arm A: two-tetrad chair G4 (MMP9-DNA-30 topology) ===")
    a = generate_arm(rng, 2, [(6, 7), (6, 7), (6, 7)])
    if a:
        cands.append(("A_chair_G4", a[1],
                      "Two-tetrad chair, the topology of MMP9-DNA-30's core -- "
                      "the only unmodified-DNA MMP-9 binder reused unchanged by "
                      "three independent groups. Novel with respect to TBA."))
        print(f"  {a[1]}  G4Hunter {a[0]:+.2f}")

    print("\n=== arm B: parallel (GGGT)n G4 (LVMH patent class) ===")
    # Three tetrads give 12 G; loops of 6-7 put the core at 30-33 nt.
    b = generate_arm(rng, 3, [(6, 7), (6, 7), (6, 7)])
    if b:
        cands.append(("B_parallel_G4", b[1],
                      "Three-tetrad parallel propeller with short loops -- the "
                      "class the LVMH MMP-9 inhibitor patent claims, and the "
                      "class our leave-one-out-validated predictor recovered "
                      "T30695 from at rank 1 of 355,914."))
        print(f"  {b[1]}  G4Hunter {b[0]:+.2f}")

    print("\n=== arm C: machine learning (Apta-MCTS vs full-length MMP-9) ===")
    if args.mcts and Path(args.mcts).exists():
        rows = list(csv.DictReader(Path(args.mcts).read_text().splitlines()))
        rows.sort(key=lambda r: -float(r["aptamer_protein_interaction_score"]))
        top = rows[0]
        # Apta-MCTS emits RNA; the user's selection is DNA, so fold U to T.
        ml_seq = top["primary_sequence"].upper().replace("U", "T")
        cands.append(("C_ml_mcts", ml_seq,
                      "Apta-MCTS, pretrained RF scorer, full-length MMP-9 "
                      "(707 aa) as target and VEGF as counter-target. INCLUDED "
                      "WITH A WARNING: we measured this scorer at AUROC 0.4927 "
                      "and 0.5055 on real aptamer-protein pairs, i.e. chance, "
                      "so this arm is here for the comparison the request asked "
                      "for, not because it is trusted."))
        print(f"  {ml_seq}  "
              f"score {top['aptamer_protein_interaction_score']} (chance-level scorer)")
    else:
        print("  MCTS output not supplied — arm C omitted")

    print("\n=== the three candidates ===\n")
    out_rows = []
    for name, seq, why in cands:
        s = score_candidate(seq)
        s["arm"] = name
        s["rationale"] = why
        out_rows.append(s)
        print(f"{name}")
        print(f"  5'-{seq}-3'   ({len(seq)} nt)")
        print(f"  G4Hunter {s['g4hunter']:+.2f}  2-tetrad motifs "
              f"{s['motif_2tetrad']}  competing MFE {s['competing_mfe']}")
        print(f"  vs TBA: {s['tba_identity']:.0%} identity, p={s['tba_p']} "
              f"| vs MMP9-DNA-30: {s['vs_MMP9_DNA_30_identity']:.0%}\n")

    payload = {
        "candidates": out_rows,
        "target": "human MMP-9, full length (UniProt P14780, 707 aa) -- matched "
                  "to SAE0077, the HEK293-expressed full-length protein used in "
                  "the user's SELEX",
        "no_affinity_predicted": (
            "No candidate carries a predicted Kd or affinity score. Four "
            "measured failures (co-folding ipTM and pose, sequence-ML "
            "protein-grouped AUROC 0.210, Apta-MCTS scorers at 0.4927/0.5055) "
            "mean any affinity number we printed would be fabricated."),
        "grading_rule": {
            "declared_before_seeing_the_selex_data": True,
            "hit": "a candidate aligns to any enriched SELEX sequence at >=70% "
                   "identity over >=20 nt",
            "strong_hit": ">=80% identity over >=30 nt to a top-10 enriched "
                          "sequence",
            "topology_hit": "the SELEX hits are G4-forming (G4Hunter > 0 with a "
                            "quadruplex motif) at a rate above a dinucleotide-"
                            "shuffled null -- this scores the PRIOR even if no "
                            "individual sequence matches",
            "null_expectation": "for 40-mers over 4^40 ~ 10^24, exact or near-"
                                "exact sequence agreement by chance is nil; the "
                                "topology_hit is the outcome that carries real "
                                "information",
            "what_would_falsify_the_prior": "SELEX hits with no G4 enrichment "
                                            "over the shuffled null. In that "
                                            "case docs/07's G4 prior and the "
                                            "focused library built on it are "
                                            "withdrawn.",
        },
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2))
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
