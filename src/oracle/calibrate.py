"""
Oracle calibration — the gate before any design run.

Before trusting the Boltz-2 co-folding oracle to drive generation, we must show
it can SEPARATE known MMP9 binders from negatives. Co-fold confidence is a
strong-but-noisy aptamer screener (ACS Synth. Biol. 2025), so this is not
optional: if separation fails, the loop would optimize noise.

Positives  : F3B (2'-F RNA aptamer, Kd ~20 nM, MMP9-selective) — the literature
             benchmark. Run as DNA too (sequence-level control).
Negatives  : random sequences of matched length + composition, and a
             poly-A control. These should score clearly lower.

Reported metrics:
  * mean/spread of oracle score for positives vs negatives
  * AUROC (rank separation)
  * a suggested decision threshold (max Youden J)

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


def build_oracle(kind: str, receptor: str | None):
    if kind == "proxy":
        return StructureProxyOracle(), "structure_proxy (CPU dry-run — NOT a binding test)"
    from oracle.boltz2 import Boltz2Oracle, Boltz2Config
    if not receptor:
        raise SystemExit("--receptor is required for the boltz2 oracle")
    cfg = Boltz2Config(receptor_sequence=read_fasta(receptor))
    return Boltz2Oracle(cfg), "boltz2_interface"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--oracle", choices=["boltz2", "proxy"], default="boltz2")
    ap.add_argument("--receptor", help="FASTA from fetch_receptor.py")
    ap.add_argument("--n-negatives", type=int, default=12)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="results/calibration.json")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    oracle, kind = build_oracle(args.oracle, args.receptor)

    # --- positives: known MMP9 binders with a published sequence -------------
    positives = []
    for apt in BENCHMARK_APTAMERS:
        if not apt.get("sequence"):
            continue
        seq_dna = apt["sequence"].replace("U", "T")   # co-fold as DNA
        positives.append((f"{apt['name']}(as DNA)", seq_dna))
    if not positives:
        raise SystemExit("no benchmark aptamer sequences available")

    # --- negatives: matched-length shuffles, randoms, poly-A ----------------
    L = len(positives[0][1])
    negatives = [(f"shuffle{i}", shuffled(positives[0][1], rng))
                 for i in range(args.n_negatives // 2)]
    negatives += [(f"random{i}", random_dna(L, rng))
                  for i in range(args.n_negatives - len(negatives) - 1)]
    negatives.append(("polyA", "A" * L))

    rows = []
    print(f"Oracle: {kind}\nScoring {len(positives)} positives and "
          f"{len(negatives)} negatives (length {L}) ...\n")
    for label, seq, cls in ([(a, b, "positive") for a, b in positives] +
                            [(a, b, "negative") for a, b in negatives]):
        sc = oracle.score(seq)
        rows.append({"label": label, "class": cls, "sequence": seq,
                     "score": sc.value, "details": sc.details})
        print(f"  {cls:8s} {label:16s} {sc.value:.4f}")

    pos = [r["score"] for r in rows if r["class"] == "positive"]
    neg = [r["score"] for r in rows if r["class"] == "negative"]
    thr, j = best_threshold(pos, neg)
    au = auroc(pos, neg)
    verdict = ("PASS — oracle separates known binders from negatives"
               if au >= 0.75 else
               "FAIL — insufficient separation; do NOT drive generation with this oracle")

    summary = {
        "oracle": kind,
        "n_positives": len(pos), "n_negatives": len(neg),
        "positive_mean": round(sum(pos) / len(pos), 4),
        "negative_mean": round(sum(neg) / len(neg), 4),
        "positive_min": round(min(pos), 4), "negative_max": round(max(neg), 4),
        "auroc": round(au, 3),
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
