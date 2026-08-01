"""
Generate the receptor MSA once and cache it.

The receptor is fixed for the whole project, but Boltz re-queries the MSA server
on every prediction unless you hand it a precomputed alignment. With hundreds to
thousands of oracle calls in the closed loop, that dominates runtime. This script
runs ONE prediction with --use_msa_server, extracts the alignment Boltz built for
the protein chain, and caches it at data/mmp9/mmp9_receptor.<ext>.

Every later call passes it via Boltz2Config(precomputed_msa=...), so no candidate
pays the MSA cost again.

Usage (on a GPU allocation, once):
    python src/target/make_msa.py --receptor data/mmp9/mmp9_receptor.fasta
"""
from __future__ import annotations
import argparse
import shutil
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from oracle.boltz2 import Boltz2Config, _write_yaml, msa_cache_base  # noqa: E402

# Boltz accepts .a3m or .csv alignments; it writes them under the run directory.
MSA_PATTERNS = ("*.a3m", "msa/*.csv", "**/msa/*.csv", "**/*.a3m")



def read_fasta(path: str) -> str:
    return "".join(l.strip() for l in Path(path).read_text().splitlines()
                   if l and not l.startswith(">"))


def find_msa(root: Path):
    """Locate the alignment Boltz generated for the protein chain."""
    seen = []
    for pat in MSA_PATTERNS:
        for hit in sorted(root.glob(pat)):
            if hit.is_file() and hit.stat().st_size > 0:
                seen.append(hit)
    # prefer the largest file — the protein MSA, not a stub for the DNA chain
    return max(seen, key=lambda p: p.stat().st_size) if seen else None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--receptor", default="data/mmp9/mmp9_receptor.fasta")
    ap.add_argument("--cache", default=None,
                    help="output path WITHOUT extension (default: keyed to the "
                         "receptor sequence so a domain switch cannot reuse the "
                         "wrong alignment)")
    ap.add_argument("--probe-dna", default="ACGTACGTACGT",
                    help="short throwaway DNA used for the one seeding run")
    ap.add_argument("--keep-workdir", action="store_true")
    args = ap.parse_args()

    receptor = read_fasta(args.receptor)
    print(f"receptor: {len(receptor)} aa from {args.receptor}")

    work = Path(tempfile.mkdtemp(prefix="msa_seed_"))
    cfg = Boltz2Config(receptor_sequence=receptor, out_root=str(work),
                       use_msa_server=True)
    yaml_path = work / "input.yaml"
    _write_yaml(yaml_path, cfg, args.probe_dna)

    import subprocess
    cmd = [cfg.boltz_bin, "predict", str(yaml_path), "--out_dir", str(work),
           "--use_msa_server", "--devices", "1", "--output_format", "pdb"]
    print("running one seeding prediction to build the MSA ...\n"
          "  (this queries the MSA server and can take several minutes)",
          flush=True)
    print("  " + " ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        print(proc.stdout[-3000:])
        raise SystemExit(f"boltz failed:\n{proc.stderr[-3000:]}")

    msa = find_msa(work)
    if msa is None:
        raise SystemExit(
            f"Could not find a generated MSA under {work}.\n"
            "Inspect that directory and pass the alignment manually via\n"
            "Boltz2Config(precomputed_msa=...). Falling back to --use_msa_server\n"
            "still works, just slower.")

    base = args.cache or str(msa_cache_base(receptor))
    dest = Path(f"{base}{msa.suffix}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(msa, dest)
    size_kb = dest.stat().st_size / 1024
    print(f"\ncached MSA: {dest}  ({size_kb:.0f} KB, from {msa.name})")
    print("Later runs pick this up automatically — no more MSA-server calls.")

    if not args.keep_workdir:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()
