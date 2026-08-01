"""
Positive control: does the oracle work on a CANONICAL protein-DNA complex?

Both readouts failed on MMP9 + aptamer (docs/06): the scalar ipTM was
noise-dominated, and the predicted pose for the known binder had zero seed-to-seed
overlap. That leaves one question unanswered, and it is the first one a reviewer
will ask:

    Is co-folding unable to screen aptamers, or is our setup broken?

This script settles it by running the SAME test on a complex the model should get
right: Zif268 (EGR1 three-finger) bound to its cognate DNA site, PDB 1AAY — a
1.6 A classic that is certainly in the training data.

  high within-seed Jaccard  -> pipeline and model are fine; the MMP9 failure is
                               specific to the aptamer / unmapped-epitope case
  low  within-seed Jaccard  -> our setup or the model is at fault, and the MMP9
                               result cannot be interpreted yet

Sequences are fetched from RCSB at runtime rather than transcribed here, because a
single wrong residue in a hand-copied control would invalidate the comparison.
The DNA is supplied as BOTH strands: Zif268 reads a duplex, and testing it against
a single strand would be a different (and unfair) question.

Usage (GPU node):
    python src/oracle/control_canonical.py --seeds 5 --out results/control_1aay.json
"""
from __future__ import annotations
import argparse
import json
import statistics as st
import subprocess
import sys
import tempfile
import time
import urllib.request
from itertools import combinations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from oracle.boltz2 import kernels_available  # noqa: E402
from pipeline.consensus import contact_residues, jaccard  # noqa: E402

RCSB_FASTA = "https://www.rcsb.org/fasta/entry/{pdb}"
PDBE_FASTA = "https://www.ebi.ac.uk/pdbe/entry/pdb/{lower}/fasta"
DNA_ALPHABET = set("ACGT")


def fetch_entities(pdb_id: str):
    """Return (protein_seqs, dna_seqs) parsed from the entry's FASTA.

    Chains are classified by alphabet rather than by header text, which varies
    between mirrors.
    """
    last = None
    for tmpl in (RCSB_FASTA, PDBE_FASTA):
        url = tmpl.format(pdb=pdb_id.upper(), lower=pdb_id.lower())
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "findaptamer/1.0"})
            with urllib.request.urlopen(req, timeout=60) as r:
                raw = r.read().decode()
        except Exception as e:
            last = e
            continue
        seqs, cur = [], []
        for line in raw.splitlines():
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                cur = []
            elif line.strip():
                cur.append(line.strip())
        if cur:
            seqs.append("".join(cur))
        prot = [s for s in seqs if not set(s.upper()) <= DNA_ALPHABET]
        dna = [s for s in seqs if set(s.upper()) <= DNA_ALPHABET and len(s) >= 4]
        if prot and dna:
            return prot, dna
        last = RuntimeError(f"parsed {len(prot)} protein / {len(dna)} DNA from {url}")
    raise SystemExit(
        f"Could not fetch {pdb_id} sequences ({last}).\n"
        "Pass them explicitly:  --protein <seq> --dna <strand1> --dna <strand2>")


def write_yaml(path: Path, protein: str, dna_strands, msa: str | None = None):
    lines = ["version: 1", "sequences:",
             "  - protein:", "      id: A", f"      sequence: {protein}"]
    if msa:
        lines.append(f"      msa: {msa}")
    for i, s in enumerate(dna_strands):
        lines += [f"  - dna:", f"      id: {chr(ord('B') + i)}",
                  f"      sequence: {s.upper()}"]
    path.write_text("\n".join(lines) + "\n")


def run_seed(protein, dna_strands, seed: int, workdir: Path):
    workdir.mkdir(parents=True, exist_ok=True)
    yml = workdir / f"seed{seed}.yaml"
    write_yaml(yml, protein, dna_strands)
    cmd = ["boltz", "predict", str(yml), "--out_dir", str(workdir),
           "--devices", "1", "--output_format", "pdb", "--override",
           "--use_msa_server", "--seed", str(seed)]
    if not kernels_available():
        cmd.append("--no_kernels")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"boltz failed (seed {seed}):\n{proc.stderr[-2000:]}")
    pdbs = sorted(workdir.rglob(f"*seed{seed}*.pdb")) or sorted(workdir.rglob("*.pdb"))
    if not pdbs:
        raise RuntimeError(f"no structure written for seed {seed} under {workdir}")
    conf = sorted(workdir.rglob(f"confidence_*seed{seed}*.json"))
    score = None
    if conf:
        d = json.loads(conf[0].read_text())
        score = d.get("iptm", d.get("complex_iptm"))
    return str(pdbs[0]), score


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pdb", default="1AAY", help="canonical protein-DNA control")
    ap.add_argument("--protein", help="override the protein sequence")
    ap.add_argument("--dna", action="append", help="override a DNA strand (repeatable)")
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--workdir", default="results/control_work")
    ap.add_argument("--out", default="results/control_1aay.json")
    args = ap.parse_args()

    if args.protein and args.dna:
        protein, dna = args.protein, args.dna
        source = "command line"
    else:
        prot_list, dna = fetch_entities(args.pdb)
        protein = prot_list[0]
        source = f"RCSB/PDBe entry {args.pdb.upper()}"
    print(f"control: {args.pdb.upper()} ({source})")
    print(f"  protein {len(protein)} aa")
    for i, s in enumerate(dna):
        print(f"  dna {chr(ord('B')+i)}: {len(s)} nt  {s}")
    print(f"  {args.seeds} seeds\n", flush=True)

    work = Path(args.workdir)
    per_seed = []
    for sd in range(args.seeds):
        t = time.time()
        struct, score = run_seed(protein, dna, sd, work / f"seed{sd}")
        cset, n_prot, n_dna = contact_residues(struct)
        per_seed.append({"seed": sd, "score": score, "structure": struct,
                         "contacts": sorted(cset)})
        print(f"  seed {sd}: contacts={len(cset)}  iptm={score}  "
              f"({time.time()-t:.0f}s)", flush=True)

    sets = [set(r["contacts"]) for r in per_seed if r["contacts"]]
    wj = (st.median([jaccard(a, b) for a, b in combinations(sets, 2)])
          if len(sets) > 1 else float("nan"))
    scores = [r["score"] for r in per_seed if r["score"] is not None]

    # MMP9 numbers to compare against (docs/06)
    MMP9 = {"positive_within_jaccard": 0.0, "negative_within_jaccard": 0.48,
            "between_sequence_jaccard": 0.351}

    if wj != wj:
        verdict = "NO DATA — contact sets could not be computed"
    elif wj >= 0.5:
        verdict = ("PIPELINE OK — the model places a canonical protein-DNA complex "
                   "reproducibly, so the MMP9 aptamer failure is specific to that "
                   "case, not an artefact of this setup")
    elif wj >= 0.3:
        verdict = ("MARGINAL — reproducibility on a canonical complex is only "
                   "moderate; interpret the MMP9 result cautiously")
    else:
        verdict = ("SETUP SUSPECT — the model cannot even place a canonical "
                   "protein-DNA complex reproducibly. Re-check the setup before "
                   "drawing conclusions about aptamers")

    summary = {
        "control": args.pdb.upper(), "seeds": args.seeds,
        "within_seed_jaccard": round(wj, 3) if wj == wj else None,
        "iptm_mean": round(st.fmean(scores), 4) if scores else None,
        "iptm_sd": round(st.pstdev(scores), 4) if len(scores) > 1 else None,
        "mmp9_comparison": MMP9,
        "verdict": verdict,
    }
    print("\n=== summary ===")
    print(json.dumps(summary, indent=2))

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"summary": summary, "per_seed": per_seed,
                               "protein": protein, "dna": dna}, indent=2))
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
