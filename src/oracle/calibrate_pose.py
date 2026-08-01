"""
Pose-reproducibility calibration — can the predicted BINDING SITE discriminate
where the scalar ipTM cannot?

Motivation (see docs/06_calibration_result_ko.md): Boltz-2's interface ipTM on
MMP9 + a 36-nt ligand has a seed-to-seed sigma of ~0.10 while the difference
between a known 20 nM binder and random sequences is ~0.003 — a signal-to-noise
ratio of 0.03, so the score is unusable. But ipTM is a single global summary; the
predicted *pose* may be far more stable than the number attached to it.

This script tests that directly, with the control that matters:

  within-sequence Jaccard   Re-predict ONE sequence with N seeds and measure the
                            median pairwise Jaccard of its receptor contact sets.
                            High = the model puts this ligand on the same patch
                            every time.

  between-sequence Jaccard  Compare the consensus patches of DIFFERENT sequences.
                            If every sequence — binder and random alike — lands on
                            the same sticky surface, within-sequence agreement is
                            high for trivial reasons and carries no specificity.

The readout is therefore not "is within-Jaccard high" but
**"is within-Jaccard higher than between-Jaccard, and is it higher for the known
binder than for the controls?"** Both must hold for pose reproducibility to be a
usable oracle signal.

Usage (GPU node):
    python src/oracle/calibrate_pose.py --receptor data/mmp9/mmp9_receptor.fasta \
        --seeds 5 --n-negatives 2 --out results/calibration_pose.json

Cost: (1 + n_negatives) x seeds predictions. Default 3 x 5 = 15 (~20 min).
"""
from __future__ import annotations
import argparse
import json
import random
import statistics as st
import sys
import time
from itertools import combinations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from target.mmp9 import BENCHMARK_APTAMERS  # noqa: E402
from pipeline.consensus import contact_residues, jaccard, residue_map, epitope_stats  # noqa: E402


def read_fasta(path: str) -> str:
    return "".join(l.strip() for l in Path(path).read_text().splitlines()
                   if l and not l.startswith(">"))


def shuffled(seq: str, rng: random.Random) -> str:
    c = list(seq)
    rng.shuffle(c)
    return "".join(c)


def predict_seeds(seq: str, chem: str, receptor_seq: str, seeds, workdir: Path):
    """Run one sequence over several seeds; return per-seed (score, contact set)."""
    from oracle.boltz2 import Boltz2Oracle, Boltz2Config
    out = []
    for sd in seeds:
        wd = workdir / f"seed{sd}"
        wd.mkdir(parents=True, exist_ok=True)
        cfg = Boltz2Config(receptor_sequence=receptor_seq, ligand_type=chem,
                           out_root=str(wd), keep_outputs=True,
                           extra_args=("--seed", str(sd)))
        t = time.time()
        sc = Boltz2Oracle(cfg).score(seq)
        struct = sc.details.get("structure")
        cset = set()
        if struct and Path(struct).exists():
            try:
                cset, _, _ = contact_residues(struct)
            except Exception as e:
                print(f"      ! contact parse failed: {e}", flush=True)
        out.append({"seed": sd, "score": sc.value, "structure": struct,
                    "contacts": sorted(cset)})
        print(f"      seed {sd}: score={sc.value:.4f}  contacts={len(cset)}  "
              f"({time.time()-t:.0f}s)", flush=True)
    return out


def consensus_patch(per_seed, min_frac=0.5):
    """Residues contacted in at least `min_frac` of the seeds."""
    sets = [set(r["contacts"]) for r in per_seed if r["contacts"]]
    if not sets:
        return set()
    need = max(1, int(len(sets) * min_frac + 0.999))
    allres = set().union(*sets)
    return {r for r in allres if sum(1 for s in sets if r in s) >= need}


def within_jaccard(per_seed):
    sets = [set(r["contacts"]) for r in per_seed if r["contacts"]]
    if len(sets) < 2:
        return float("nan")
    return st.median([jaccard(a, b) for a, b in combinations(sets, 2)])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--receptor", required=True)
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--n-negatives", type=int, default=2)
    ap.add_argument("--seed", type=int, default=0, help="RNG seed for the controls")
    ap.add_argument("--workdir", default="results/pose_work")
    ap.add_argument("--out", default="results/calibration_pose.json")
    args = ap.parse_args()

    receptor_seq = read_fasta(args.receptor)
    resmap = residue_map(args.receptor)
    if len(resmap) != len(receptor_seq):
        raise SystemExit(
            f"receptor length {len(receptor_seq)} != residue map {len(resmap)}; "
            f"re-run src/target/fetch_receptor.py to embed the range spec.")

    rng = random.Random(args.seed)
    seeds = list(range(args.seeds))
    workroot = Path(args.workdir)
    workroot.mkdir(parents=True, exist_ok=True)

    # positive: the one MMP9 binder with a published sequence, in its real chemistry
    apt = next(a for a in BENCHMARK_APTAMERS if a.get("sequence"))
    chem = "rna" if "RNA" in apt["type"].upper() else "dna"
    subjects = [(apt["name"], apt["sequence"], "positive")]
    for i in range(args.n_negatives):
        subjects.append((f"shuffle{i}", shuffled(apt["sequence"], rng), "negative"))
    subjects.append(("polyA", "A" * len(apt["sequence"]), "negative"))

    print(f"receptor {len(receptor_seq)} aa (UniProt {resmap[0]}-{resmap[-1]}), "
          f"chemistry {chem}")
    print(f"{len(subjects)} sequences x {len(seeds)} seeds = "
          f"{len(subjects)*len(seeds)} predictions\n")

    records = []
    for name, seq, cls in subjects:
        print(f"  {cls:8s} {name}", flush=True)
        per_seed = predict_seeds(seq, chem, receptor_seq, seeds, workroot / name)
        patch = consensus_patch(per_seed)
        wj = within_jaccard(per_seed)
        scores = [r["score"] for r in per_seed]
        est = epitope_stats(patch, resmap, None)
        rec = {"name": name, "class": cls, "sequence": seq, "chemistry": chem,
               "score_mean": round(st.fmean(scores), 4),
               "score_sd": round(st.pstdev(scores), 4) if len(scores) > 1 else 0.0,
               "within_jaccard": round(wj, 3) if wj == wj else None,
               "consensus_patch_uniprot": sorted(resmap[c - 1] for c in patch
                                                 if 1 <= c <= len(resmap)),
               **est, "per_seed": per_seed}
        records.append(rec)
        print(f"    -> within-seed Jaccard {rec['within_jaccard']}, "
              f"consensus patch {len(patch)} residues, "
              f"score {rec['score_mean']:.4f}+/-{rec['score_sd']:.4f}\n", flush=True)

    # ---- the control: do DIFFERENT sequences land on the SAME patch? ---------
    patches = {r["name"]: set(r["consensus_patch_uniprot"]) for r in records}
    between = []
    for a, b in combinations(records, 2):
        j = jaccard(patches[a["name"]], patches[b["name"]])
        between.append({"pair": [a["name"], b["name"]], "jaccard": round(j, 3)})

    pos = [r["within_jaccard"] for r in records
           if r["class"] == "positive" and r["within_jaccard"] is not None]
    neg = [r["within_jaccard"] for r in records
           if r["class"] == "negative" and r["within_jaccard"] is not None]
    bet = [b["jaccard"] for b in between]

    within_all = [r["within_jaccard"] for r in records
                  if r["within_jaccard"] is not None]
    specificity = (st.fmean(within_all) - st.fmean(bet)) if within_all and bet else float("nan")

    verdict = []
    if not within_all:
        verdict.append("NO DATA — no contact sets could be computed")
    else:
        if specificity == specificity and specificity < 0.1:
            verdict.append(
                "NOT SPECIFIC — different sequences land on the same patch as often "
                "as one sequence lands on its own; pose agreement is a property of "
                "the receptor surface, not of the ligand")
        if pos and neg and st.fmean(pos) <= st.fmean(neg):
            verdict.append(
                "NOT DISCRIMINATIVE — the known binder is no more reproducible "
                "than the controls")
        if not verdict:
            verdict.append(
                "PROMISING — the binder's pose is both self-consistent and more "
                "specific than chance; pose reproducibility may serve as an oracle")

    summary = {
        "n_sequences": len(records), "seeds": len(seeds),
        "positive_within_jaccard": round(st.fmean(pos), 3) if pos else None,
        "negative_within_jaccard": round(st.fmean(neg), 3) if neg else None,
        "between_sequence_jaccard": round(st.fmean(bet), 3) if bet else None,
        "specificity_margin": round(specificity, 3) if specificity == specificity else None,
        "note": ("within-sequence Jaccard must exceed between-sequence Jaccard "
                 "(specificity) AND be higher for the positive (discrimination). "
                 "Either alone is not enough."),
        "verdict": " | ".join(verdict),
    }

    print("\n=== summary ===")
    print(json.dumps(summary, indent=2))
    print("\nbetween-sequence patch overlap:")
    for b in between:
        print(f"  {b['pair'][0]:12s} vs {b['pair'][1]:12s}  J={b['jaccard']}")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"summary": summary, "between": between,
                               "records": records}, indent=2, default=str))
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
