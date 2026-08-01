"""
Does a sequence-based aptamer-protein model actually READ the protein?

This is the gate before any in-silico SELEX. The plan is to throw random 40-mers
at MMP-9 and let a model select binders across rounds. That plan only means
something if the model's score depends on WHICH protein it is scoring against.
If it does not, then "binds MMP-9" is not a prediction about MMP-9 at all -- the
same sequences would win against any target, and the SELEX simulation would be
an elaborate way of ranking DNA by base composition.

There is a specific reason to suspect exactly that. Published aptamer-protein
interaction sets build their NEGATIVES by randomly re-pairing aptamers and
proteins. Under that construction a model can score near-perfectly by learning
"is this string a real aptamer" and never looking at the protein column at all:
real aptamers are positives, and a randomly paired real aptamer is a negative
only by virtue of the pairing, which the model is free to ignore... except it
cannot distinguish those two cases from the aptamer alone. So the honest
question is empirical, and it is answered by ablation.

Four conditions, same model, same folds:
  full            aptamer features + protein features
  aptamer_only    protein features removed
  protein_only    aptamer features removed
  protein_shuffled  protein features permuted across rows, destroying the
                    pairing while preserving every marginal distribution

The decisive comparison is full vs aptamer_only, and full vs protein_shuffled.
If removing or scrambling the protein does not hurt performance, the model is
protein-blind and unusable as a target-specific oracle.

LEAKAGE CONTROL: the same aptamer appears in many pairs. A random split puts
copies of one aptamer on both sides and inflates every number. Folds here are
grouped by aptamer identity, so an aptamer seen in training never appears in
test. The ungrouped number is also reported, because the gap between them is
itself worth seeing.

Data: github.com/nedaemami/AptaNet (Sci Rep 2021), features pre-extracted by
the authors -- 340 aptamer k-mer frequencies, 300 protein PseAAC descriptors.

Usage:
    python src/ml/ablation.py --data <path>/AptaNet/Dataset.csv
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold, StratifiedKFold


def load(path: str):
    d = pd.read_csv(path)
    y = d["Class"].to_numpy()
    apt_cols = [c for c in d.columns if c.startswith("aptamer_")]
    pro_cols = [c for c in d.columns if c != "Class" and c not in apt_cols]
    return d, y, apt_cols, pro_cols


def aptamer_groups(d: pd.DataFrame, apt_cols: list[str]) -> np.ndarray:
    """Group id per distinct aptamer, so folds never split one across sides.

    Aptamer identity is recovered from its feature vector: the k-mer profile is
    a deterministic function of the sequence, so identical rows are the same
    aptamer even though the sequence itself is not in the file.
    """
    key = [hash(t) for t in map(tuple, d[apt_cols].round(6).to_numpy())]
    uniq = {k: i for i, k in enumerate(dict.fromkeys(key))}
    return np.array([uniq[k] for k in key])


def evaluate(X, y, groups, seed=0, n_splits=5, grouped=True):
    """Mean out-of-fold AUROC and accuracy."""
    if grouped:
        splitter = GroupKFold(n_splits=n_splits).split(X, y, groups)
    else:
        splitter = StratifiedKFold(n_splits=n_splits, shuffle=True,
                                   random_state=seed).split(X, y)
    aucs, accs = [], []
    for tr, te in splitter:
        if len(np.unique(y[te])) < 2:
            continue
        clf = RandomForestClassifier(n_estimators=300, n_jobs=-1,
                                     random_state=seed, class_weight="balanced")
        clf.fit(X[tr], y[tr])
        p = clf.predict_proba(X[te])[:, 1]
        aucs.append(roc_auc_score(y[te], p))
        accs.append((clf.predict(X[te]) == y[te]).mean())
    return float(np.mean(aucs)), float(np.mean(accs))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data", required=True, help="AptaNet Dataset.csv")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="results/ml_ablation.json")
    args = ap.parse_args()

    d, y, apt_cols, pro_cols = load(args.data)
    groups = aptamer_groups(d, apt_cols)
    rng = np.random.default_rng(args.seed)

    print(f"{len(d)} pairs | {y.sum()} positive / {(y == 0).sum()} negative")
    print(f"{len(apt_cols)} aptamer features, {len(pro_cols)} protein features")
    print(f"{len(np.unique(groups))} distinct aptamers "
          f"({len(d) / len(np.unique(groups)):.1f} pairs each)\n")

    Xa = d[apt_cols].to_numpy()
    Xp = d[pro_cols].to_numpy()
    # Permute protein rows: pairing destroyed, marginals identical.
    Xp_shuf = Xp[rng.permutation(len(Xp))]

    conditions = {
        "full": np.hstack([Xa, Xp]),
        "aptamer_only": Xa,
        "protein_only": Xp,
        "protein_shuffled": np.hstack([Xa, Xp_shuf]),
    }

    rows = {}
    print(f"{'condition':<20} {'AUROC (grouped)':>16} {'acc':>7} "
          f"{'AUROC (random split)':>21}")
    print("-" * 70)
    for name, X in conditions.items():
        auc_g, acc_g = evaluate(X, y, groups, args.seed, grouped=True)
        auc_r, _ = evaluate(X, y, groups, args.seed, grouped=False)
        rows[name] = {"auroc_grouped": round(auc_g, 4),
                      "accuracy_grouped": round(acc_g, 4),
                      "auroc_random_split": round(auc_r, 4)}
        print(f"{name:<20} {auc_g:>16.4f} {acc_g:>7.4f} {auc_r:>21.4f}")

    full = rows["full"]["auroc_grouped"]
    apt = rows["aptamer_only"]["auroc_grouped"]
    shuf = rows["protein_shuffled"]["auroc_grouped"]
    drop_removed = full - apt
    drop_shuffled = full - shuf

    # A model that reads the protein must lose ground when the protein is taken
    # away or scrambled. 0.02 AUROC is a deliberately lenient bar.
    reads_protein = drop_removed > 0.02 and drop_shuffled > 0.02
    verdict = (
        "USABLE as a target-specific oracle: removing the protein costs "
        f"{drop_removed:+.4f} AUROC and scrambling it costs {drop_shuffled:+.4f}, "
        "so the score does depend on which protein is being scored."
        if reads_protein else
        "PROTEIN-BLIND -- NOT usable for MMP-9. Removing the protein costs "
        f"{drop_removed:+.4f} AUROC and scrambling it costs {drop_shuffled:+.4f}. "
        "The model scores the aptamer alone, so it would return the same "
        "ranking for MMP-9 as for any other target. An in-silico SELEX driven "
        "by this score would not be selecting for MMP-9 binding."
    )

    print(f"\nprotein removed : {drop_removed:+.4f} AUROC")
    print(f"protein shuffled: {drop_shuffled:+.4f} AUROC")
    print(f"\nVERDICT: {verdict}")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({
        "n_pairs": int(len(d)), "n_positive": int(y.sum()),
        "n_distinct_aptamers": int(len(np.unique(groups))),
        "conditions": rows,
        "delta_auroc_protein_removed": round(drop_removed, 4),
        "delta_auroc_protein_shuffled": round(drop_shuffled, 4),
        "reads_protein": bool(reads_protein),
        "verdict": verdict,
        "source": "github.com/nedaemami/AptaNet (Sci Rep 2021)",
    }, indent=2))
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
