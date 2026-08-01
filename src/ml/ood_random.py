"""
Can the model tell a RANDOM 40-mer from a real binder? (out-of-distribution gate)

The in-silico SELEX plan is: generate random 40-nt DNA, score it against MMP-9,
keep the winners, mutate, repeat. For that to select anything, the scorer must
rank real binders above random DNA.

There is a structural reason to doubt it, and it is not the protein-blindness
question that `ablation.py` already settled. It is the negative class. In the
AptaNet-style dataset every negative is a REAL APTAMER paired with the wrong
protein. The model has therefore never once been shown a random sequence. Asked
"does this random 40-mer bind?", it is extrapolating off the edge of its
training distribution, and a model can be perfectly good at its training task
while being useless at that one.

This script measures it directly. For each protein in the data that has known
binders, it scores:

    positives   the aptamers that really bind this protein
    negatives   the real aptamers paired with it as dataset negatives
    random40    freshly generated random 40-mers  <- never seen in training

and reports where random DNA lands. The number that matters is
AUROC(positives vs random40). If random 40-mers score like real binders, a
SELEX loop driven by this model would be selecting noise, and we should not
run it.

Training uses aptamer-grouped folds so no test aptamer was seen in training.

Usage:
    python src/ml/ood_random.py --data <path>/AptaNet/Dataset.csv
"""
from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold

KMER_K = (1, 2, 3, 4)          # 4 + 16 + 64 + 256 = 340, matching the dataset
BASES = "ACGT"


def kmer_features(seq: str) -> np.ndarray:
    """AptaNet's aptamer encoding: normalised k-mer frequencies for k=1..4.

    Column order must match the dataset exactly, so it is generated the same
    way the column names are: alphabetically within each k, k ascending.
    """
    from itertools import product
    seq = seq.upper().replace("U", "T")
    out = []
    for k in KMER_K:
        counts = {"".join(p): 0 for p in product(BASES, repeat=k)}
        n = max(len(seq) - k + 1, 1)
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i + k]
            if kmer in counts:
                counts[kmer] += 1
        out.extend(counts[m] / n for m in sorted(counts))
    return np.array(out)


def random_dna(n: int, rng: random.Random) -> str:
    return "".join(rng.choice(BASES) for _ in range(n))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data", required=True)
    ap.add_argument("--length", type=int, default=40,
                    help="length of the random library sequences")
    ap.add_argument("--n-random", type=int, default=200)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="results/ml_ood_random.json")
    args = ap.parse_args()

    d = pd.read_csv(args.data)
    y = d["Class"].to_numpy()
    apt_cols = [c for c in d.columns if c.startswith("aptamer_")]
    pro_cols = [c for c in d.columns if c != "Class" and c not in apt_cols]
    Xa, Xp = d[apt_cols].to_numpy(), d[pro_cols].to_numpy()
    X = np.hstack([Xa, Xp])

    # Verify our featuriser reproduces the dataset's column semantics before
    # trusting any score built on it: a mismatch in column order would be
    # invisible and would silently corrupt every result below.
    assert len(kmer_features("ACGT")) == len(apt_cols), (
        f"featuriser emits {len(kmer_features('ACGT'))} columns, "
        f"dataset has {len(apt_cols)}")

    groups_key = [hash(t) for t in map(tuple, d[apt_cols].round(6).to_numpy())]
    uniq = {k: i for i, k in enumerate(dict.fromkeys(groups_key))}
    groups = np.array([uniq[k] for k in groups_key])

    prot_key = [hash(t) for t in map(tuple, d[pro_cols].round(6).to_numpy())]
    puniq = {k: i for i, k in enumerate(dict.fromkeys(prot_key))}
    prot_id = np.array([puniq[k] for k in prot_key])
    print(f"{len(d)} pairs | {len(np.unique(groups))} aptamers "
          f"| {len(np.unique(prot_id))} proteins\n")

    rng = random.Random(args.seed)
    rand_seqs = [random_dna(args.length, rng) for _ in range(args.n_random)]
    Xr_apt = np.vstack([kmer_features(s) for s in rand_seqs])

    # Train on grouped folds; score the held-out test aptamers AND the random
    # library with the same fold model, so nothing the random set is compared
    # against was trained on.
    rows, aurocs, rand_scores_all, pos_scores_all = [], [], [], []
    for tr, te in GroupKFold(n_splits=5).split(X, y, groups):
        clf = RandomForestClassifier(n_estimators=300, n_jobs=-1,
                                     random_state=args.seed,
                                     class_weight="balanced")
        clf.fit(X[tr], y[tr])

        for pid in np.unique(prot_id[te]):
            sel = te[(prot_id[te] == pid)]
            pos = sel[y[sel] == 1]
            if len(pos) < 2:
                continue
            pvec = Xp[sel[0]]
            s_pos = clf.predict_proba(np.hstack(
                [Xa[pos], np.tile(pvec, (len(pos), 1))]))[:, 1]
            s_rand = clf.predict_proba(np.hstack(
                [Xr_apt, np.tile(pvec, (len(Xr_apt), 1))]))[:, 1]
            lab = np.r_[np.ones(len(s_pos)), np.zeros(len(s_rand))]
            auc = roc_auc_score(lab, np.r_[s_pos, s_rand])
            aurocs.append(auc)
            pos_scores_all.extend(s_pos.tolist())
            rand_scores_all.extend(s_rand.tolist())
            rows.append({"protein": int(pid), "n_binders": int(len(pos)),
                         "mean_binder_score": round(float(s_pos.mean()), 4),
                         "mean_random_score": round(float(s_rand.mean()), 4),
                         "auroc_binder_vs_random": round(float(auc), 4)})

    mean_auc = float(np.mean(aurocs)) if aurocs else float("nan")
    mp, mr = float(np.mean(pos_scores_all)), float(np.mean(rand_scores_all))
    print(f"{'protein':>8} {'binders':>8} {'binder':>8} {'random':>8} {'AUROC':>7}")
    print("-" * 45)
    for r in sorted(rows, key=lambda r: -r["auroc_binder_vs_random"])[:15]:
        print(f"{r['protein']:>8} {r['n_binders']:>8} "
              f"{r['mean_binder_score']:>8.3f} {r['mean_random_score']:>8.3f} "
              f"{r['auroc_binder_vs_random']:>7.3f}")

    print(f"\n{len(rows)} proteins tested")
    print(f"mean AUROC (real binder vs random {args.length}-mer): {mean_auc:.4f}")
    print(f"mean score  binders {mp:.4f}   random {mr:.4f}   "
          f"gap {mp - mr:+.4f}")

    usable = mean_auc > 0.65
    verdict = (
        f"USABLE for in-silico SELEX: real binders outrank random "
        f"{args.length}-mers at AUROC {mean_auc:.3f}, so selection pressure is real."
        if usable else
        f"NOT USABLE for in-silico SELEX: real binders separate from random "
        f"{args.length}-mers at only AUROC {mean_auc:.3f}. Every negative this "
        f"model was trained on was itself a real aptamer, so it was never "
        f"taught what random DNA looks like. Rounds driven by this score would "
        f"select whatever base composition the model happens to favour, not "
        f"MMP-9 binding."
    )
    print(f"\nVERDICT: {verdict}")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({
        "random_length": args.length, "n_random": args.n_random,
        "n_proteins_tested": len(rows),
        "mean_auroc_binder_vs_random": round(mean_auc, 4),
        "mean_binder_score": round(mp, 4),
        "mean_random_score": round(mr, 4),
        "per_protein": rows, "usable": bool(usable), "verdict": verdict,
    }, indent=2))
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
