"""
SELEX library and primer design for MMP-9: 5'-[FWD 20]-[N40]-[revcomp REV 20]-3'

Standard primer design rules (Tm, GC, hairpin, dimer, 3'-end) are necessary
here but not sufficient, because this selection has a constraint a generic
primer tool does not know about.

THE PROJECT-SPECIFIC CONSTRAINT: the aptamer class we are selecting for is
G-quadruplex. All three MMP-9 aptamers with a reported fold are G4 or G-quartet
(docs/07), and the counter-example F3B scores -1.17 by G4Hunter. If the fixed
primer regions are G-rich they will (a) fold into their own G4, and (b) supply
G-tracts that combine with G-runs in the random region to form quadruplexes
that span the primer/N40 boundary. Either way the constant regions -- present
in every single member of the library -- start competing for the target, and a
constant region that binds is the fastest way to ruin a SELEX: it enriches
uniformly and carries no sequence information.

So primers here are additionally required to be G4-dead: no GG run at all, a
G4Hunter score at or below zero, and no 2-tetrad motif anywhere in the full
construct including across both junctions.

The second project-specific check is the junction, and it has to be phrased
carefully. Asking "does the assembled construct contain a G4 motif" is vacuous:
the stand-in random region is poly-G and a 40-nt G run contains 30 motifs on its
own, so that test fires for every pair. The question that matters is whether a
motif RECRUITS constant-region bases -- a quadruplex inside the random region is
the whole point, one borrowing a G from a primer is a liability present in every
member. motif_uses_constant_region() locates each motif and rejects only spans
that cross a boundary.

Everything else is the standard battery:
  strand    the G4 tests run on the sequence that lands IN the library, which
            for the reverse primer is revcomp(REV), not REV
  Tm        nearest-neighbour, SantaLucia 1998 unified parameters, salt- and
            concentration-corrected; the pair must match within --tm-tol
  GC        40-60%
  3' end    no more than 2 G/C in the last 5 (a strong 3' clamp promotes
            mispriming), and no 3'-terminal T
  hairpin   ViennaRNA under DNA parameters
  dimer     self- and cross-dimer, with the 3' end weighted because only a
            3'-anchored dimer is extendable by polymerase
  repeats   no run of 4+ identical bases, no long dinucleotide repeat

Usage:
    python src/design/selex_library.py                 # design and report
    python src/design/selex_library.py --n-library 20  # + example members
"""
from __future__ import annotations

import argparse
import json
import math
import random
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from analysis.g4 import g4_motifs, g4hunter_windows  # noqa: E402

COMP = str.maketrans("ACGT", "TGCA")

# SantaLucia (1998) unified nearest-neighbour parameters.
# dH kcal/mol, dS cal/(mol.K)
NN = {
    "AA": (-7.9, -22.2), "AT": (-7.2, -20.4), "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7), "GT": (-8.4, -22.4), "CT": (-7.8, -21.0),
    "GA": (-8.2, -22.2), "CG": (-10.6, -27.2), "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9),
}
INIT_GC = (0.1, -2.8)      # initiation with terminal G/C
INIT_AT = (2.3, 4.1)       # initiation with terminal A/T
R = 1.987                  # cal/(mol.K)


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def _nn(pair: str):
    if pair in NN:
        return NN[pair]
    return NN[revcomp(pair)]


def tm(seq: str, conc_nM: float = 250.0, na_mM: float = 50.0,
       mg_mM: float = 1.5) -> float:
    """Nearest-neighbour Tm with a monovalent-equivalent salt correction.

    Mg2+ is folded into an equivalent Na+ using the usual 120*sqrt(Mg) rule;
    it matters because PCR buffers are Mg-dominated and ignoring it
    underestimates Tm by several degrees.
    """
    seq = seq.upper()
    if len(seq) < 2:
        return float("nan")
    dh, ds = 0.0, 0.0
    for i in range(len(seq) - 1):
        a, b = _nn(seq[i:i + 2])
        dh += a
        ds += b
    for end in (seq[0], seq[-1]):
        a, b = INIT_GC if end in "GC" else INIT_AT
        dh += a
        ds += b
    na_eq = (na_mM + 120.0 * (max(mg_mM, 0.0) ** 0.5)) / 1000.0
    ds += 0.368 * (len(seq) - 1) * math.log(na_eq)
    ct = conc_nM * 1e-9
    return (dh * 1000.0) / (ds + R * math.log(ct / 4.0)) - 273.15


def gc(seq: str) -> float:
    return (seq.count("G") + seq.count("C")) / len(seq)


def max_run(seq: str) -> int:
    best = run = 1
    for i in range(1, len(seq)):
        run = run + 1 if seq[i] == seq[i - 1] else 1
        best = max(best, run)
    return best


def has_dinuc_repeat(seq: str, n: int = 4) -> bool:
    """(XY)n repeats slip during PCR and are a common SELEX artefact."""
    for i in range(len(seq) - 2 * n + 1):
        unit = seq[i:i + 2]
        if seq[i:i + 2 * n] == unit * n:
            return True
    return False


_DNA_PARAMS_LOADED = False


def hairpin_dG(seq: str):
    """MFE under DNA parameters. The parameter table is loaded once -- doing it
    per call, which is the obvious way to write this, dominated the runtime
    since primer sampling calls this hundreds of thousands of times."""
    global _DNA_PARAMS_LOADED
    try:
        import RNA
    except ImportError:
        return None
    if not _DNA_PARAMS_LOADED:
        RNA.params_load_DNA_Mathews2004()
        _DNA_PARAMS_LOADED = True
    _, mfe = RNA.fold_compound(seq, RNA.md()).mfe()
    return round(mfe, 2)


def dimer_score(a: str, b: str, three_prime_window: int = 6) -> int:
    """Longest complementary overlap, counted only when it anchors a 3' end.

    A dimer that does not involve a 3' terminus cannot be extended by
    polymerase, so it is far less damaging and is not scored here.
    """
    rb = b.translate(COMP)[::-1]
    best = 0
    for off in range(-(len(a) - 1), len(rb)):
        match = span = 0
        for i in range(len(a)):
            j = i - off
            if 0 <= j < len(rb):
                span += 1
                if a[i] == rb[j]:
                    match += 1
        if span >= 4 and match == span:
            tail_a = off + len(rb) >= len(a) - three_prime_window
            if tail_a:
                best = max(best, span)
    return best


def g4_dead(seq: str) -> tuple[bool, dict]:
    """Primer must not form, or contribute to, a quadruplex.

    The junction test uses poly-G flanks as the worst case: N40 is random, so
    some members WILL present G-tracts next to the primer, and we need the
    primer to be inert even then.
    """
    score, _, _ = g4hunter_windows(seq)
    m = g4_motifs(seq)
    flanked = "G" * 6 + seq + "G" * 6
    mj = g4_motifs(flanked)
    ok = (max_run_g(seq) < 2 and score <= 0.0
          and not m["two_tetrad"] and not m["three_tetrad"]
          and not mj["two_tetrad"] and not mj["three_tetrad"])
    return ok, {"g4hunter": round(score, 3),
                "max_G_run": max_run_g(seq),
                "motif_alone": len(m["two_tetrad"]),
                "motif_at_junction": len(mj["two_tetrad"])}


def motif_uses_constant_region(construct: str, len5: int, len_n: int) -> bool:
    """Does a quadruplex motif in the assembled library use CONSTANT bases?

    The obvious version of this check -- build the worst case as
    FWD + G*40 + revcomp(REV) and ask whether a motif exists -- is vacuous: a
    40-nt poly-G run contains 30 motifs entirely by itself, so it fires for
    every primer pair and rejects all of them. That is a property of the
    stand-in random region, not of the primers.

    What actually matters is whether a motif RECRUITS constant-region bases. A
    quadruplex living wholly inside the random region is the point of the
    library; one that borrows a G from a primer is a constant-region liability,
    because that G is present in every single member. The realistic junction
    risk is small but real: a primer ending in a lone G, next to a random region
    starting with G, forms a GG that straddles the boundary.

    So find each motif's span and reject only spans that cross a boundary.
    """
    import re
    lo, hi = len5, len5 + len_n          # half-open span of the random region
    for m in re.finditer(r"(?=(G{2,}\w{1,7}G{2,}\w{1,7}G{2,}\w{1,7}G{2,}))",
                         construct):
        start = m.start()
        end = start + len(m.group(1))
        if start < lo or end > hi:
            return True
    return False


def max_run_g(seq: str) -> int:
    best = run = 0
    for c in seq:
        run = run + 1 if c == "G" else 0
        best = max(best, run)
    return best


def screen_primer(p: str, tm_lo: float, tm_hi: float, hp_floor: float,
                  library_strand: str | None = None):
    """Screen a primer. `library_strand` is the sequence that actually ends up
    IN the library member, which is NOT always the primer.

    This distinction is the whole point. The aptamer strand is
    FWD + N40 + revcomp(REV), so for the reverse primer it is revcomp(REV) that
    sits in every library member. Screening REV itself for quadruplexes lets the
    exact failure we are guarding against in through the back door: G4-dead
    requires no GG run, which together with GC 40-60% pushes a primer C-rich,
    and the reverse complement of a C-rich primer is G-rich. So the G4 test is
    applied to the library strand while the PCR tests (3' end, self-dimer,
    hairpin) stay on the primer, because those are properties of the oligo the
    polymerase actually uses.
    """
    in_library = library_strand if library_strand is not None else p
    t = tm(p)
    if not (tm_lo <= t <= tm_hi):
        return None, f"Tm {t:.1f} outside [{tm_lo},{tm_hi}]"
    if not (0.40 <= gc(p) <= 0.60):
        return None, f"GC {gc(p):.0%}"
    if max_run(p) > 3:
        return None, "homopolymer run > 3"
    if has_dinuc_repeat(p):
        return None, "dinucleotide repeat"
    if p[-1] == "T":
        return None, "3'-terminal T"
    if sum(1 for c in p[-5:] if c in "GC") > 2:
        return None, "3' GC clamp too strong"
    ok, g4 = g4_dead(in_library)
    if not ok:
        return None, "G4-forming or G4-contributing"
    hp = hairpin_dG(p)
    if hp is not None and hp < hp_floor:
        return None, f"hairpin {hp}"
    if dimer_score(p, p) >= 5:
        return None, "self-dimer at 3'"
    return {"sequence": p, "in_library": in_library, "tm": round(t, 1),
            "gc": round(gc(p), 3), "hairpin_dG": hp, **g4}, "pass"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--random-len", type=int, default=40)
    ap.add_argument("--primer-len", type=int, default=20)
    ap.add_argument("--tm-lo", type=float, default=58.0)
    ap.add_argument("--tm-hi", type=float, default=64.0)
    ap.add_argument("--tm-tol", type=float, default=1.5,
                    help="max Tm difference between the two primers")
    ap.add_argument("--hairpin-floor", type=float, default=-2.0,
                    help="reject primers with a hairpin below this (kcal/mol)")
    ap.add_argument("--tries", type=int, default=400000)
    ap.add_argument("--n-library", type=int, default=0,
                    help="also emit this many example library members")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="results/selex_design.json")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    # A/T-leaning sampling: the G4-dead requirement rejects G-rich primers, so
    # sampling uniformly wastes most draws.
    alphabet = "AATTCCG"
    fwd_c, rev_c, rejects = [], [], {}
    for _ in range(args.tries):
        p = "".join(rng.choice(alphabet) for _ in range(args.primer_len))
        # forward: the primer IS the 5' constant region of the library
        rec, why = screen_primer(p, args.tm_lo, args.tm_hi, args.hairpin_floor)
        if rec:
            fwd_c.append(rec)
        else:
            rejects[why.split()[0]] = rejects.get(why.split()[0], 0) + 1
        # reverse: revcomp(primer) is what lands in the library
        rec2, why2 = screen_primer(p, args.tm_lo, args.tm_hi,
                                   args.hairpin_floor,
                                   library_strand=revcomp(p))
        if rec2:
            rev_c.append(rec2)
        if len(fwd_c) >= 250 and len(rev_c) >= 250:
            break

    print(f"{len(fwd_c)} forward / {len(rev_c)} reverse primers passed")
    for k, v in sorted(rejects.items(), key=lambda x: -x[1])[:6]:
        print(f"  rejected {v:>7}  {k}")
    if not fwd_c or not rev_c:
        raise SystemExit("not enough candidates; loosen --tm-lo/--tm-hi")

    # Pair: matched Tm, no cross-dimer, and the assembled construct must stay
    # G4-clean across both junctions.
    best = None
    for f in fwd_c:
        for r in rev_c:
            if f["sequence"] == r["sequence"]:
                continue
            if abs(f["tm"] - r["tm"]) > args.tm_tol:
                continue
            if dimer_score(f["sequence"], r["sequence"]) >= 5:
                continue
            construct = f["in_library"] + "G" * args.random_len + r["in_library"]
            if motif_uses_constant_region(construct, len(f["in_library"]),
                                          args.random_len):
                continue
            # Objective: matched Tm first, then the least G-rich primers as
            # OLIGOS. The library strand is already guaranteed G4-dead, but
            # making the 3' constant region G-poor forces its reverse
            # complement -- the REV primer -- to be G-rich, and that is
            # unavoidable rather than a bug. A primer with four G-tracts could
            # form an intramolecular quadruplex and prime badly; with three it
            # cannot, though intermolecular association is still possible at
            # high concentration. So among Tm-matched pairs, prefer the one
            # whose primers carry the fewest and shortest G-tracts.
            def g_liability(seq):
                runs = [len(x) for x in re.findall(r"G{2,}", seq)]
                return len(runs) + 0.25 * sum(runs)
            pen = (abs(f["tm"] - r["tm"])
                   + 0.5 * (g_liability(f["sequence"])
                            + g_liability(r["sequence"])))
            if best is None or pen < best[0]:
                best = (pen, f, r)
    if best is None:
        raise SystemExit("no compatible primer pair found")

    _, fwd, rev = best
    template5 = fwd["sequence"]
    template3 = revcomp(rev["sequence"])
    total = args.primer_len * 2 + args.random_len

    print(f"\n=== SELEX library ({total} nt) ===\n")
    print(f"  5'-{template5}-(N){args.random_len}-{template3}-3'\n")
    print(f"  FWD primer  5'-{fwd['sequence']}-3'   "
          f"Tm {fwd['tm']}  GC {fwd['gc']:.0%}  G4Hunter {fwd['g4hunter']}")
    print(f"  REV primer  5'-{rev['sequence']}-3'   "
          f"Tm {rev['tm']}  GC {rev['gc']:.0%}  G4Hunter {rev['g4hunter']}")
    print(f"\n  Tm difference {abs(fwd['tm']-rev['tm']):.1f} C   "
          f"cross-dimer {dimer_score(fwd['sequence'], rev['sequence'])} bp   "
          f"max G-run  fwd {fwd['max_G_run']} / rev {rev['max_G_run']}")
    print(f"  diversity of the random region: 4^{args.random_len} = "
          f"10^{args.random_len * 0.602:.0f}")

    payload = {
        "construct": f"5'-{template5}-(N){args.random_len}-{template3}-3'",
        "total_length": total, "random_length": args.random_len,
        "fwd_primer": fwd, "rev_primer": rev,
        "constant_5prime": template5, "constant_3prime": template3,
        "tm_difference": round(abs(fwd["tm"] - rev["tm"]), 2),
        "cross_dimer_bp": dimer_score(fwd["sequence"], rev["sequence"]),
        "design_notes": {
            "g4_dead_primers": (
                "Primers are required to have no GG run, G4Hunter <= 0, and no "
                "quadruplex motif even when flanked by poly-G. The MMP-9 "
                "aptamers with a known fold are all G4 (docs/07), so a G-rich "
                "constant region would compete for the target in every library "
                "member and enrich without carrying sequence information."),
            "tm_model": ("SantaLucia 1998 unified nearest-neighbour, 250 nM "
                         "primer, 50 mM Na+, 1.5 mM Mg2+ as monovalent "
                         "equivalent. Re-check against your polymerase's "
                         "buffer before ordering."),
            "ssDNA_generation": (
                "PCR gives duplex; SELEX needs the sense (aptamer) strand "
                "alone. Recommended: order the REVERSE primer 5'-phosphorylated "
                "and treat the PCR product with lambda exonuclease, which "
                "digests the 5'-phosphorylated strand -- that is the antisense "
                "strand, since it is the one primed by REV -- and leaves "
                "FWD+N40+revcomp(REV) intact. This avoids the streptavidin "
                "carryover and NaOH-denaturation losses of the biotin route. "
                "The forward primer needs no modification. If you prefer "
                "biotin-streptavidin instead, biotinylate REV, not FWD."),
            "counter_selection": (
                "SAE0077 is HEK293-expressed full-length MMP-9, so the FnII "
                "gelatin-binding exosite is present. ssDNA binds FnII and "
                "ACTIVATES MMP-9 (Shimada 2018), and MMP-2 carries the same "
                "FnII domains -- so an FnII binder would likely cross-react "
                "with MMP-2 and be useless as an MMP-9-specific diagnostic. "
                "Counter-select against MMP-2 and against gelatin/collagen."),
        },
    }

    if args.n_library:
        lib = []
        for i in range(args.n_library):
            n = "".join(rng.choice("ACGT") for _ in range(args.random_len))
            lib.append({"id": f"N{args.random_len}_{i:04d}",
                        "sequence": template5 + n + template3,
                        "random_region": n})
        payload["example_members"] = lib
        print(f"\n  {args.n_library} example members in the JSON")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2))
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
