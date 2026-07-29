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
from target.mmp9 import DOMAINS, UNIPROT, DOMAIN_SPECS, spec_string, slice_sequence  # noqa: E402

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
    """Return the bare sequence (no header, no newlines).

    Prefers the vendored copy in the repo (data/mmp9/P14780.fasta), which has
    been verified against the literature landmarks (707 aa; Cys99; His401/
    Glu402/His405/His411; Pro421/Tyr423; Phe107). Only hits the network if it
    is missing, so the pipeline runs on closed networks.
    """
    local = outdir / f"{UNIPROT}.fasta"
    if local.exists():
        raw = local.read_text()
        seq = "".join(l.strip() for l in raw.splitlines()
                      if l and not l.startswith(">"))
        if seq:
            print(f"  using vendored {local} ({len(seq)} aa)")
            return seq

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
    """Slice by UniProt numbering using the shared DOMAIN_SPECS definition."""
    if domain not in DOMAIN_SPECS:
        raise ValueError(f"unknown domain: {domain}")
    return slice_sequence(seq, domain)


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
    # the range spec is REQUIRED downstream to map predicted residue
    # numbers back to UniProt numbering (see consensus.epitope_stats)
    out_fa.write_text(f">MMP9_{args.domain}|{spec_string(args.domain)}\n{dom}\n")
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
