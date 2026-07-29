"""
Stage 4 — consensus validation of the top candidates.

A single co-fold confidence score is a weak aptamer screener. This stage asks a
harder question of each surviving candidate:

    "Do independent predictions agree on WHERE and HOW WELL it binds?"

For every candidate we re-predict the complex several times — multiple seeds of
Boltz-2, and (if installed) Chai-1 as an independently trained second model —
then score three things:

  score_mean / score_std   confidence and its stability across predictions
  contact_jaccard          agreement of the receptor-residue contact set between
                           predictions (median pairwise Jaccard). This is the
                           real test: a candidate that lands on a different patch
                           each time is a false positive no matter how confident
                           any single run looks.
  epitope_fraction         fraction of contacts falling on the intended epitope
                           (and, for MMP9, how much lands on the FnII exosite —
                           reported so an inhibitor programme can counter-select)

Candidates are ranked by a consensus score that requires BOTH confidence and
agreement, so a lucky single high score cannot carry a candidate through.

Usage (GPU node):
  python src/pipeline/consensus.py \
      --receptor data/mmp9/mmp9_receptor.fasta \
      --candidates results/pilot/top_candidates.csv \
      --top 20 --seeds 3 --use-chai --out results/consensus
"""
from __future__ import annotations
import argparse
import csv
import json
import statistics
import sys
from itertools import combinations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from target.mmp9 import DOMAINS, EPITOPES  # noqa: E402

# residue-level contact cutoff between an aptamer atom and a receptor atom
CONTACT_CUTOFF_A = 5.0


# ---------------------------------------------------------------------------
# structure parsing (PDB or mmCIF) — minimal, dependency-light
# ---------------------------------------------------------------------------
def _parse_structure(path: Path):
    """Return {chain_id: [(resseq, resname, [(x,y,z), ...]), ...]}."""
    from Bio.PDB import PDBParser, MMCIFParser
    parser = MMCIFParser(QUIET=True) if path.suffix.lower() in (".cif", ".mmcif") \
        else PDBParser(QUIET=True)
    st = parser.get_structure("m", str(path))
    model = next(iter(st))
    out = {}
    for chain in model:
        residues = []
        for res in chain:
            coords = [tuple(a.get_coord()) for a in res]
            if coords:
                residues.append((res.id[1], res.get_resname().strip(), coords))
        if residues:
            out[chain.id] = residues
    return out


def _is_nucleic(resnames) -> bool:
    nuc = {"DA", "DT", "DG", "DC", "A", "U", "G", "C", "T",
           "DA3", "DA5", "DT3", "DT5", "DG3", "DG5", "DC3", "DC5"}
    hits = sum(1 for r in resnames if r in nuc)
    return hits >= max(1, len(resnames) // 2)


def contact_residues(structure_path: str, cutoff: float = CONTACT_CUTOFF_A):
    """Receptor residue numbers contacted by the nucleic-acid chain.

    Returns (contact_set, n_protein_res, n_dna_res). Empty set if the file does
    not contain both a protein and a nucleic-acid chain.
    """
    chains = _parse_structure(Path(structure_path))
    prot, dna = {}, {}
    for cid, residues in chains.items():
        if _is_nucleic([rn for _, rn, _ in residues]):
            dna[cid] = residues
        else:
            prot[cid] = residues
    if not prot or not dna:
        return set(), 0, 0

    import numpy as np
    dna_atoms = np.array([c for residues in dna.values() for _, _, coords in residues
                          for c in coords], dtype=float)
    cut2 = float(cutoff) ** 2
    contacts = set()
    for residues in prot.values():
        for resseq, _, coords in residues:
            P = np.asarray(coords, dtype=float)
            d2 = ((P[:, None, :] - dna_atoms[None, :, :]) ** 2).sum(-1)
            if bool((d2 <= cut2).any()):
                contacts.add(resseq)
    n_prot = sum(len(r) for r in prot.values())
    n_dna = sum(len(r) for r in dna.values())
    return contacts, n_prot, n_dna


def jaccard(a: set, b: set) -> float:
    if not a and not b:
        return 0.0
    return len(a & b) / len(a | b)


# ---------------------------------------------------------------------------
# receptor numbering: predictions are 1..N over the sliced receptor, but the
# epitope definitions use UniProt numbering. Map with the slice offset.
# ---------------------------------------------------------------------------
def slice_offset(receptor_fasta: str) -> int:
    """First UniProt residue number of the receptor slice (from its header)."""
    header = Path(receptor_fasta).read_text().splitlines()[0]
    # e.g. ">MMP9_catalytic_107-443"
    for tok in header.replace(">", "").split("_"):
        if "-" in tok:
            a, b = tok.split("-", 1)
            if a.isdigit() and b.isdigit():
                return int(a)
    return DOMAINS["catalytic_domain"][0]


def epitope_stats(contacts_local: set, offset: int, epitope_key: str | None):
    """Translate local numbering to UniProt and measure epitope / FnII overlap."""
    uni = {c + offset - 1 for c in contacts_local}
    fa, fb = DOMAINS["fnII_inserts"]
    fn_hits = sum(1 for r in uni if fa <= r <= fb)
    out = {"n_contacts": len(uni),
           "fnII_fraction": round(fn_hits / len(uni), 3) if uni else 0.0}
    if epitope_key and epitope_key in EPITOPES:
        target = set(EPITOPES[epitope_key]["residues"])
        near = sum(1 for r in uni if any(abs(r - t) <= 8 for t in target))
        out["epitope"] = epitope_key
        out["epitope_fraction"] = round(near / len(uni), 3) if uni else 0.0
    return out


# ---------------------------------------------------------------------------
def read_candidates(path: str, top: int):
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f):
            rows.append(r)
    rows = rows[:top]
    return [(r.get("sequence"), float(r.get("oracle_score", 0) or 0)) for r in rows]


def build_predictors(receptor_seq: str, seeds: int, use_chai: bool, outroot: Path):
    """Return [(name, callable(seq, workdir) -> (score, structure_path_or_None))]."""
    preds = []
    from oracle.boltz2 import Boltz2Oracle, Boltz2Config

    for s in range(seeds):
        def make(seed=s):
            def run(seq, workdir):
                cfg = Boltz2Config(receptor_sequence=receptor_seq,
                                   out_root=str(workdir),
                                   extra_args=("--seed", str(seed)))
                sc = Boltz2Oracle(cfg).score(seq)
                struct = None
                for pat in ("*.pdb", "*.cif"):
                    hits = sorted(Path(workdir).rglob(pat))
                    if hits:
                        struct = str(hits[0])
                        break
                return sc.value, struct, sc.details
            return run
        preds.append((f"boltz2_seed{s}", make()))

    if use_chai:
        from oracle.chai1 import Chai1Oracle, Chai1Config

        def run_chai(seq, workdir):
            cfg = Chai1Config(receptor_sequence=receptor_seq,
                              out_root=str(workdir), seed=0)
            sc = Chai1Oracle(cfg).score(seq)
            return sc.value, sc.details.get("structure"), sc.details
        preds.append(("chai1", run_chai))
    return preds


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--receptor", required=True)
    ap.add_argument("--candidates", default="results/pilot/top_candidates.csv")
    ap.add_argument("--top", type=int, default=20)
    ap.add_argument("--seeds", type=int, default=3, help="Boltz-2 seeds per candidate")
    ap.add_argument("--use-chai", action="store_true",
                    help="add Chai-1 as an independent second model")
    ap.add_argument("--epitope", default=None,
                    choices=[None, *EPITOPES.keys()],
                    help="report contact overlap with this epitope")
    ap.add_argument("--min-score", type=float, default=0.5)
    ap.add_argument("--min-jaccard", type=float, default=0.5)
    ap.add_argument("--out", default="results/consensus")
    args = ap.parse_args()

    receptor_seq = "".join(l.strip() for l in Path(args.receptor).read_text().splitlines()
                           if l and not l.startswith(">"))
    offset = slice_offset(args.receptor)
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)

    cands = read_candidates(args.candidates, args.top)
    predictors = build_predictors(receptor_seq, args.seeds, args.use_chai, outdir)
    print(f"{len(cands)} candidates x {len(predictors)} predictions "
          f"(receptor {len(receptor_seq)} aa, UniProt offset {offset})\n")

    results = []
    for i, (seq, pilot_score) in enumerate(cands, 1):
        scores, contact_sets, per_pred = [], [], {}
        for name, run in predictors:
            wd = outdir / f"cand{i:03d}" / name
            wd.mkdir(parents=True, exist_ok=True)
            try:
                val, struct, details = run(seq, wd)
            except Exception as e:
                print(f"  ! {name} failed on candidate {i}: {e}")
                continue
            scores.append(val)
            cset = set()
            if struct and Path(struct).exists():
                try:
                    cset, _, _ = contact_residues(struct)
                except Exception as e:
                    print(f"  ! contact parse failed ({name}): {e}")
            if cset:
                contact_sets.append(cset)
            per_pred[name] = {"score": val, "structure": struct,
                              "n_contacts": len(cset)}

        if not scores:
            print(f"{i:3d}. {seq[:28]}...  ALL PREDICTIONS FAILED")
            continue

        mean = statistics.fmean(scores)
        std = statistics.pstdev(scores) if len(scores) > 1 else 0.0
        jac = (statistics.median([jaccard(a, b)
                                  for a, b in combinations(contact_sets, 2)])
               if len(contact_sets) > 1 else 0.0)
        # union of contacts seen in a majority of predictions
        stable = set()
        if contact_sets:
            need = max(1, len(contact_sets) // 2 + 1)
            allres = set().union(*contact_sets)
            stable = {r for r in allres
                      if sum(1 for cs in contact_sets if r in cs) >= need}
        est = epitope_stats(stable, offset, args.epitope)

        # consensus score: confidence AND agreement must both hold
        consensus = mean * (0.5 + 0.5 * jac)
        passed = (mean >= args.min_score and jac >= args.min_jaccard)

        rec = {"rank_in": i, "sequence": seq, "length": len(seq),
               "pilot_score": pilot_score,
               "score_mean": round(mean, 4), "score_std": round(std, 4),
               "n_predictions": len(scores),
               "contact_jaccard": round(jac, 3),
               "consensus_score": round(consensus, 4),
               "stable_contacts_uniprot": sorted(r + offset - 1 for r in stable),
               **est, "pass": passed, "per_predictor": per_pred}
        results.append(rec)
        flag = "PASS" if passed else "    "
        print(f"{i:3d}. {seq[:26]}..  mean={mean:.3f}±{std:.3f}  "
              f"jaccard={jac:.2f}  consensus={consensus:.3f}  {flag}")

    results.sort(key=lambda r: r["consensus_score"], reverse=True)
    (outdir / "consensus.json").write_text(json.dumps(results, indent=2, default=str))
    with open(outdir / "consensus_ranked.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["rank", "sequence", "length", "consensus_score", "score_mean",
                    "score_std", "contact_jaccard", "n_contacts", "fnII_fraction",
                    "pass"])
        for k, r in enumerate(results, 1):
            w.writerow([k, r["sequence"], r["length"], r["consensus_score"],
                        r["score_mean"], r["score_std"], r["contact_jaccard"],
                        r.get("n_contacts", 0), r.get("fnII_fraction", 0),
                        r["pass"]])

    n_pass = sum(1 for r in results if r["pass"])
    print(f"\n{n_pass}/{len(results)} passed "
          f"(score>={args.min_score}, jaccard>={args.min_jaccard})")
    print(f"-> {outdir}/consensus_ranked.csv")
    if n_pass == 0 and results:
        print("\nNo candidate reached consensus. Options: relax --min-jaccard, "
              "increase --seeds, or widen the pilot search (more rounds).")


if __name__ == "__main__":
    main()
