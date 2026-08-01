"""
Oracle calibration — the gate before any design run.

Before trusting the Boltz-2 co-folding oracle to drive generation, we must show
it can SEPARATE known MMP9 binders from negatives. Co-fold confidence is a
strong-but-noisy aptamer screener (ACS Synth. Biol. 2025), so this is not
optional: if separation fails, the loop would optimize noise.

Positives  : F3B (Kd ~20 nM, MMP9-selective). It is a 2'-F pyrimidine RNA
             aptamer, so it is co-folded AS RNA — folding it as B-DNA would be a
             different molecule and a correct oracle could rank it below random.
Negatives  : shuffles of the positive, random sequences of matched length, and a
             poly-A control — all in the positive's chemistry, so the test
             measures binding rather than chemistry.

Reported metrics:
  * per-sequence mean over --seeds predictions, plus the seed-to-seed spread
  * AUROC (rank separation) and a separation z-score vs the negative distribution
  * a suggested decision threshold (max Youden J)

IMPORTANT LIMITATION: only one MMP9 binder has a published sequence, so this gate
rests on a single positive. AUROC is quantized to 1/n_negatives and noisy. Treat
a marginal result as inconclusive rather than as validation.

Usage (on the GPU server, after fetch_receptor.py):
    python3 src/oracle/calibrate.py --receptor data/mmp9/mmp9_receptor.fasta \
        --n-negatives 12 --out results/calibration.json

Add --oracle proxy to dry-run the harness on CPU (no Boltz, no binding signal —
it only checks the plumbing; it will NOT validate the binding oracle).
"""
from __future__ import annotations
import argparse
import json
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from target.mmp9 import BENCHMARK_APTAMERS  # noqa: E402
from oracle.interface import StructureProxyOracle  # noqa: E402


def read_fasta(path: str) -> str:
    lines = Path(path).read_text().splitlines()
    return "".join(l.strip() for l in lines if l and not l.startswith(">"))


def shuffled(seq: str, rng: random.Random) -> str:
    chars = list(seq)
    rng.shuffle(chars)
    return "".join(chars)


def random_dna(n: int, rng: random.Random) -> str:
    return "".join(rng.choice("AGCT") for _ in range(n))


def auroc(pos, neg) -> float:
    """Rank-based AUROC; ties count as 0.5."""
    if not pos or not neg:
        return float("nan")
    wins = 0.0
    for p in pos:
        for n in neg:
            wins += 1.0 if p > n else (0.5 if p == n else 0.0)
    return wins / (len(pos) * len(neg))


def best_threshold(pos, neg):
    """Threshold maximizing Youden's J (TPR - FPR)."""
    cands = sorted(set(pos + neg))
    best, best_j = None, -1.0
    for t in cands:
        tpr = sum(1 for p in pos if p >= t) / len(pos)
        fpr = sum(1 for n in neg if n >= t) / len(neg)
        j = tpr - fpr
        if j > best_j:
            best, best_j = t, j
    return best, best_j


def build_oracle(kind: str, receptor: str | None, ligand_type: str = "dna",
                 seed: int | None = None):
    """Oracle for one chemistry (dna|rna) and optionally one seed."""
    if kind == "proxy":
        return StructureProxyOracle(), "structure_proxy (CPU dry-run — NOT a binding test)"
    from oracle.boltz2 import Boltz2Oracle, Boltz2Config
    if not receptor:
        raise SystemExit("--receptor is required for the boltz2 oracle")
    extra = ("--seed", str(seed)) if seed is not None else ()
    cfg = Boltz2Config(receptor_sequence=read_fasta(receptor),
                       ligand_type=ligand_type, extra_args=extra)
    return Boltz2Oracle(cfg), "boltz2_interface"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--oracle", choices=["boltz2", "proxy"], default="boltz2")
    ap.add_argument("--receptor", help="FASTA from fetch_receptor.py")
    ap.add_argument("--n-negatives", type=int, default=12)
    ap.add_argument("--seeds", type=int, default=3,
                    help="predictions per sequence; >1 estimates the "
                         "oracle run-to-run spread")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="results/calibration.json")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    kind = ("structure_proxy (CPU dry-run — NOT a binding test)"
            if args.oracle == "proxy" else "boltz2_interface")

    # --- positives: known MMP9 binders, scored as the chemistry they ARE ------
    # F3B is a 2'-F pyrimidine RNA aptamer; co-folding it as B-DNA is a different
    # molecule with a different (A-form) fold, so a CORRECT oracle could rank it
    # below random DNA. Declare the real entity type instead.
    positives = []
    for apt in BENCHMARK_APTAMERS:
        if not apt.get("sequence"):
            continue
        chem = "rna" if "RNA" in apt["type"].upper() else "dna"
        positives.append((apt["name"], apt["sequence"], chem))
    if not positives:
        raise SystemExit("no benchmark aptamer sequences available")

    # --- negatives: matched-length shuffles, randoms, poly-A ----------------
    # Negatives share the positive's chemistry and length so the comparison is
    # like-for-like (an RNA positive vs DNA negatives would test chemistry, not
    # binding).
    p_name, p_seq, p_chem = positives[0]
    L = len(p_seq)
    negatives = [(f"shuffle{i}", shuffled(p_seq, rng))
                 for i in range(args.n_negatives // 2)]
    negatives += [(f"random{i}", random_dna(L, rng))
                  for i in range(args.n_negatives - len(negatives) - 1)]
    negatives.append(("polyA", "A" * L))

    rows = []
    seeds = list(range(args.seeds))
    print(f"Oracle: {kind}\n{len(positives)} positive(s) x {len(seeds)} seed(s) "
          f"and {len(negatives)} negatives (length {L}, chemistry {p_chem}) ...\n")

    import time as _time
    total_calls = (len(positives) + len(negatives)) * len(seeds)
    state = {"done": 0, "t0": _time.time()}

    def score_one(label, seq, chem, cls):
        vals = []
        for sd in seeds:
            t = _time.time()
            o, _ = build_oracle(args.oracle, args.receptor, chem,
                                sd if args.oracle != "proxy" else None)
            sc = o.score(seq)
            vals.append(sc.value)
            state["done"] += 1
            # Each co-fold takes minutes and its subprocess output is captured,
            # so without a per-seed line the log looks frozen for a long time.
            el = _time.time() - state["t0"]
            eta = el / state["done"] * (total_calls - state["done"])
            print(f"    [{state['done']:>3d}/{total_calls}] {label}#s{sd} "
                  f"= {sc.value:.4f}  ({_time.time()-t:.0f}s, ETA {eta/60:.0f}m)",
                  flush=True)
            rows.append({"label": f"{label}#s{sd}", "class": cls, "chemistry": chem,
                         "sequence": seq, "seed": sd, "score": sc.value,
                         "details": sc.details})
        m = sum(vals) / len(vals)
        spread = (max(vals) - min(vals)) if len(vals) > 1 else 0.0
        print(f"  {cls:8s} {label:16s} {m:.4f}"
              + (f"  (spread {spread:.3f} over {len(vals)} seeds)" if len(vals) > 1 else ""),
              flush=True)
        return m

    pos = [score_one(n, s, c, "positive") for n, s, c in positives]
    neg = [score_one(n, s, p_chem, "negative") for n, s in negatives]
    thr, j = best_threshold(pos, neg)
    au = auroc(pos, neg)
    # With a single published positive, AUROC is quantized to 1/len(neg) and very
    # noisy. Report a separation z-score against the negative distribution too,
    # and say plainly how weak the evidence is.
    import statistics as _st
    neg_mean = _st.fmean(neg)
    neg_sd = _st.pstdev(neg) if len(neg) > 1 else 0.0
    z = ((_st.fmean(pos) - neg_mean) / neg_sd) if neg_sd > 1e-9 else float("nan")
    strong = (au >= 0.75) and (z >= 1.0 if z == z else True)
    verdict = ("PASS — oracle separates the known binder from negatives"
               if strong else
               "FAIL — insufficient separation; do NOT drive generation with this oracle")

    summary = {
        "oracle": kind,
        "n_positives": len(pos), "n_negatives": len(neg),
        "positive_mean": round(sum(pos) / len(pos), 4),
        "negative_mean": round(sum(neg) / len(neg), 4),
        "positive_min": round(min(pos), 4), "negative_max": round(max(neg), 4),
        "auroc": round(au, 3),
        "separation_z": round(z, 2) if z == z else None,
        "n_positive_sequences": len(positives),
        "seeds": len(seeds),
        "evidence_note": ("only one MMP9 binder has a published sequence, so this "
                          "gate rests on a single positive — treat a marginal "
                          "result as inconclusive, not as validation"),
        "suggested_threshold": round(thr, 4) if thr is not None else None,
        "youden_j": round(j, 3),
        "verdict": verdict,
    }
    print("\n" + json.dumps(summary, indent=2))
    if args.oracle == "proxy":
        print("\nNOTE: the proxy oracle has no protein in it. This run only "
              "verifies the harness works — it does NOT validate binding.")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"summary": summary, "rows": rows}, indent=2))
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
