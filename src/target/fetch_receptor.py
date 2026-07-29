"""
Fetch the MMP9 receptor sequence and structures.

Run this ON THE SERVER (open network). The dev sandbox blocks UniProt/RCSB, so
nothing is vendored into the repo — the pipeline fails loudly with a clear
message rather than using a transcribed (error-prone) sequence.

Downloads:
  * UniProt P14780 canonical FASTA           -> data/mmp9/P14780.fasta
  * PDB structures (default 1GKC, 5TH6, 1L6J) -> data/mmp9/<ID>.pdb
  * A trimmed catalytic-domain FASTA          -> data/mmp9/mmp9_catalytic.fasta

The catalytic-domain slice is what goes into the Boltz-2 receptor field.
Two slicing options (see src/target/mmp9.py for the domain map):
  --domain catalytic      residues 107-443 (catalytic incl. FnII inserts)
  --domain catalytic_nofn residues 107-215 + 391-443 joined (FnII removed;
                          smaller -> much cheaper co-folding, and avoids the
                          FnII surface entirely)
"""
from __future__ import annotations
import argparse
import json
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from target.mmp9 import DOMAINS, UNIPROT  # noqa: E402

UNIPROT_URLS = [
    f"https://rest.uniprot.org/uniprotkb/{UNIPROT}.fasta",
    f"https://www.uniprot.org/uniprot/{UNIPROT}.fasta",
    f"https://www.ebi.ac.uk/proteins/api/proteins/{UNIPROT}",
]
PDB_URLS = [
    "https://files.rcsb.org/download/{pdb}.pdb",
    "https://www.ebi.ac.uk/pdbe/entry-files/download/pdb{lower}.ent",
]


def _get(url: str, timeout=60) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "findaptamer/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def fetch_uniprot_fasta(outdir: Path) -> str:
    """Return the bare sequence (no header, no newlines)."""
    last = None
    for url in UNIPROT_URLS:
        try:
            raw = _get(url).decode()
        except Exception as e:                      # try the next mirror
            last = e
            continue
        if url.endswith(".fasta"):
            lines = [l for l in raw.splitlines() if l and not l.startswith(">")]
            seq = "".join(lines)
            (outdir / f"{UNIPROT}.fasta").write_text(raw)
        else:                                        # EBI proteins API (JSON)
            seq = json.loads(raw)["sequence"]["sequence"]
            (outdir / f"{UNIPROT}.fasta").write_text(f">{UNIPROT}\n{seq}\n")
        if seq:
            return seq
    raise RuntimeError(f"Could not fetch {UNIPROT} from any source. Last: {last}")


def fetch_pdb(pdb_id: str, outdir: Path):
    last = None
    for tmpl in PDB_URLS:
        url = tmpl.format(pdb=pdb_id.upper(), lower=pdb_id.lower())
        try:
            data = _get(url)
        except Exception as e:
            last = e
            continue
        out = outdir / f"{pdb_id.upper()}.pdb"
        out.write_bytes(data)
        return out
    print(f"  ! {pdb_id}: download failed ({last})", file=sys.stderr)
    return None


def slice_domain(seq: str, domain: str) -> str:
    """Slice by UniProt numbering (1-based, inclusive)."""
    if domain == "catalytic":
        a, b = DOMAINS["catalytic_domain"]
        return seq[a - 1:b]
    if domain == "catalytic_nofn":
        ca, cb = DOMAINS["catalytic_domain"]
        fa, fb = DOMAINS["fnII_inserts"]
        return seq[ca - 1:fa - 1] + seq[fb:cb]
    if domain == "full":
        return seq
    raise ValueError(f"unknown domain: {domain}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", default="data/mmp9")
    ap.add_argument("--domain", default="catalytic",
                    choices=["catalytic", "catalytic_nofn", "full"])
    ap.add_argument("--pdb", nargs="*", default=["1GKC", "5TH6", "1L6J"])
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Fetching UniProt {UNIPROT} ...")
    seq = fetch_uniprot_fasta(outdir)
    print(f"  full length: {len(seq)} aa (expected 707)")
    if len(seq) != 707:
        print("  ! WARNING: unexpected length; check the entry/isoform",
              file=sys.stderr)

    dom = slice_domain(seq, args.domain)
    out_fa = outdir / "mmp9_receptor.fasta"
    out_fa.write_text(f">MMP9_{args.domain}\n{dom}\n")
    print(f"  {args.domain}: {len(dom)} aa -> {out_fa}")

    for pdb in args.pdb:
        p = fetch_pdb(pdb, outdir)
        if p:
            print(f"  {pdb} -> {p}")

    # sanity check: catalytic zinc motif should be present in the full sequence
    motif_ok = seq[400:411].startswith("H")   # His401 (1-based) => index 400
    print(json.dumps({"uniprot_len": len(seq), "receptor_len": len(dom),
                      "receptor_fasta": str(out_fa),
                      "His401_check": motif_ok}, indent=2))


if __name__ == "__main__":
    main()
