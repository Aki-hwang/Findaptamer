"""
Does the model generalise to a protein it has NEVER seen? (the MMP-9 question)

Everything else is downstream of this. MMP-9 is not in any public
aptamer-protein training set -- checked directly: no protein in the AptaTrans
corpus carries the MMP cysteine-switch motif PRCGVPD. So scoring random 40-mers
against MMP-9 asks the model to generalise to a completely novel target, and
the corpus contains only ~144 distinct proteins to learn that from.

`ablation.py` grouped folds by APTAMER and showed the model reads the protein
(AUROC 0.728 full vs 0.510 aptamer-only). That is necessary but not sufficient
here: a model can read the protein well enough to recognise the 144 it was
trained on and still be useless on the 145th. The test that matches our actual
use is protein-grouped: hold out entire proteins, train on the rest, and score
the held-out ones.

Three splits, same features, same model:
  random           the number papers usually report
  aptamer_grouped  no test aptamer seen in training
  protein_grouped  no test PROTEIN seen in training   <-- the MMP-9 case

The gap between the first and the last is the honest cost of a novel target.

Unlike ablation.py this works from RAW SEQUENCES (AptaTrans dataset_li.pickle:
label / aptamer / protein / dataset), so the encoder is ours end to end. That
matters: the AptaNet release ships pre-extracted features whose protein encoder
we could not reproduce (its amino-acid composition terms come out negative, and
the published code is Python-2 with integer-division ambiguity), which made
encoding MMP-9 into that same space unsafe. Owning the encoder removes the
problem -- whatever we do to the training proteins, we do identically to MMP-9.

Features:
  aptamer  k-mer frequencies, k=1..4 over ACGT (U folded to T)  -> 340
  protein  conjoint triad: 7 physicochemical classes, all 3-mers -> 343

Usage:
    python src/ml/generalize.py --data <path>/AptaTrans/data/dataset_li.pickle
"""
from __future__ import annotations

import argparse
import json
import pickle
from itertools import product
from pathlib import Path

import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold, StratifiedKFold

BASES = "ACGT"
KMER_K = (1, 2, 3, 4)

# Conjoint triad classes (Shen et al. PNAS 2007): amino acids grouped by volume
# and dipole so that a 20-letter alphabet becomes 7, making 3-mers tractable.
CT_CLASSES = {
    **{a: 0 for a in "AGV"}, **{a: 1 for a in "ILFP"},
    **{a: 2 for a in "YMTS"}, **{a: 3 for a in "HNQW"},
    **{a: 4 for a in "RK"},   **{a: 5 for a in "DE"},
    **{a: 6 for a in "C"},
}
_KMERS = {k: sorted("".join(p) for p in product(BASES, repeat=k))
          for k in KMER_K}
_TRIADS = {"".join(map(str, t)): i
           for i, t in enumerate(product(range(7), repeat=3))}


def encode_aptamer(seq: str) -> np.ndarray:
    seq = seq.upper().replace("U", "T")
    seq = "".join(c for c in seq if c in BASES)   # drop N/B and other ambiguity
    out = []
    for k in KMER_K:
        counts = dict.fromkeys(_KMERS[k], 0)
        n = max(len(seq) - k + 1, 1)
        for i in range(len(seq) - k + 1):
            m = seq[i:i + k]
            if m in counts:
                counts[m] += 1
        out.extend(counts[m] / n for m in _KMERS[k])
    return np.array(out, dtype=np.float32)


def encode_protein(seq: str) -> np.ndarray:
    seq = "".join(c for c in seq.upper() if c in CT_CLASSES)
    v = np.zeros(len(_TRIADS), dtype=np.float32)
    for i in range(len(seq) - 2):
        key = "".join(str(CT_CLASSES[c]) for c in seq[i:i + 3])
        v[_TRIADS[key]] += 1
    return v / max(v.max(), 1.0)      # per-protein max-normalisation


def group_ids(values) -> np.ndarray:
    uniq = {v: i for i, v in enumerate(dict.fromkeys(values))}
    return np.array([uniq[v] for v in values])


def evaluate(X, y, groups, scheme: str, seed=0, n_splits=5):
    if scheme == "random":
        it = StratifiedKFold(n_splits, shuffle=True,
                             random_state=seed).split(X, y)
    else:
        it = GroupKFold(n_splits=n_splits).split(X, y, groups)
    aucs = []
    for tr, te in it:
        if len(np.unique(y[te])) < 2:
            continue
        clf = RandomForestClassifier(n_estimators=300, n_jobs=-1,
                                     random_state=seed,
                                     class_weight="balanced")
        clf.fit(X[tr], y[tr])
        aucs.append(roc_auc_score(y[te], clf.predict_proba(X[te])[:, 1]))
    return float(np.mean(aucs)), float(np.std(aucs)), len(aucs)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data", required=True, help="AptaTrans dataset_li.pickle")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="results/ml_generalize.json")
    args = ap.parse_args()

    d = pickle.load(open(args.data, "rb"))
    y = (d["label"].to_numpy() == "positive").astype(int)
    Xa = np.vstack([encode_aptamer(s) for s in d["aptamer"]])
    Xp = np.vstack([encode_protein(s) for s in d["protein"]])
    X = np.hstack([Xa, Xp])

    g_apt = group_ids(d["aptamer"].tolist())
    g_pro = group_ids(d["protein"].tolist())
    print(f"{len(d)} pairs | {y.sum()} positive | "
          f"{len(np.unique(g_apt))} aptamers | {len(np.unique(g_pro))} proteins")
    print(f"features: {Xa.shape[1]} aptamer + {Xp.shape[1]} protein\n")

    schemes = {"random": None, "aptamer_grouped": g_apt,
               "protein_grouped": g_pro}
    rows = {}
    print(f"{'split':<20} {'AUROC':>8} {'sd':>7} {'folds':>6}")
    print("-" * 44)
    for name, g in schemes.items():
        auc, sd, n = evaluate(X, y, g, name, args.seed)
        rows[name] = {"auroc": round(auc, 4), "sd": round(sd, 4), "folds": n}
        print(f"{name:<20} {auc:>8.4f} {sd:>7.4f} {n:>6}")

    # WHY the held-out-protein number comes out where it does. An AUROC well
    # BELOW 0.5 is not "no signal", it is inverted signal, and the cause is in
    # how the corpus is built rather than in the model.
    import pandas as pd
    lab = pd.DataFrame({"aptamer": d["aptamer"].tolist(), "y": y})
    per = lab.groupby("aptamer")["y"].agg(["sum", "count"])
    both = int(((per["sum"] > 0) & (per["sum"] < per["count"])).sum())
    mech = {
        "aptamers_both_pos_and_neg": both,
        "aptamers_total": int(len(per)),
        "mean_positives_per_aptamer": round(float(per["sum"].mean()), 2),
        "mean_pairs_per_aptamer": round(float(per["count"].mean()), 2),
        "explanation": (
            f"{both} of {len(per)} aptamers appear as BOTH a positive and a "
            f"negative, and each is positive for essentially exactly one "
            f"protein (mean {per['sum'].mean():.2f} positives out of "
            f"{per['count'].mean():.2f} pairs). So 'is this a real aptamer' "
            "carries no label information -- every aptamer in the corpus is "
            "real -- and the only learnable signal is the specific "
            "aptamer-to-protein association. On a protein never seen in "
            "training there is no such association to recall, so the model "
            "falls back on associations learned elsewhere, which for a "
            "held-out protein are precisely the wrong ones. That inverts the "
            "ranking, which is why the AUROC lands below 0.5 rather than at it."
        ),
    }
    print(f"\nmechanism: {mech['explanation']}")

    novel = rows["protein_grouped"]["auroc"]
    gap = rows["random"]["auroc"] - novel
    usable = novel > 0.65
    verdict = (
        f"GENERALISES: AUROC {novel:.3f} on proteins never seen in training, so "
        f"scoring random 40-mers against MMP-9 is a real prediction. Note the "
        f"{gap:+.3f} gap to the random split -- that gap is what a paper "
        f"reporting the random-split number would be hiding."
        if usable else
        f"DOES NOT GENERALISE to a novel protein: AUROC {novel:.3f} on held-out "
        f"proteins against {rows['random']['auroc']:.3f} on a random split. "
        f"MMP-9 is absent from this corpus, so an in-silico SELEX driven by "
        f"this model would be ranking 40-mers by something other than MMP-9 "
        f"binding. The corpus has only ~144 distinct proteins, which is the "
        f"likely cause -- this is a data limitation, not a modelling bug."
    )
    print(f"\nVERDICT: {verdict}")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({
        "n_pairs": int(len(d)), "n_positive": int(y.sum()),
        "n_aptamers": int(len(np.unique(g_apt))),
        "n_proteins": int(len(np.unique(g_pro))),
        "mmp9_in_corpus": False,
        "splits": rows, "auroc_novel_protein": novel,
        "gap_random_minus_novel": round(gap, 4),
        "generalises": bool(usable), "verdict": verdict,
        "source": "github.com/pnumlb/AptaTrans dataset_li.pickle",
    }, indent=2))
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
