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

# Boltz's YAML `msa:` field takes an .a3m or .csv alignment. Which of those it
# leaves on disk (and where) varies by version, so search the whole run tree by
# extension rather than guessing a fixed layout.
MSA_SUFFIXES = (".a3m", ".csv")


def read_fasta(path: str) -> str:
    return "".join(l.strip() for l in Path(path).read_text().splitlines()
                   if l and not l.startswith(">"))


def find_msa(root: Path, verbose: bool = True):
    """Locate the alignment Boltz generated for the protein chain.

    Returns the largest reusable alignment found, or None. On failure it lists
    what IS in the tree, because the alternative — silently falling back to the
    MSA server for every prediction — turns a 20-minute calibration into hours.
    """
    seen = [p for p in root.rglob("*")
            if p.is_file() and p.suffix.lower() in MSA_SUFFIXES
            and p.stat().st_size > 0]
    if seen:
        best = max(seen, key=lambda p: p.stat().st_size)
        if verbose:
            print(f"  candidate alignments: "
                  + ", ".join(f"{p.relative_to(root)}({p.stat().st_size}B)"
                              for p in sorted(seen)[:8]))
        return best
    if verbose:
        files = sorted(p for p in root.rglob("*") if p.is_file())
        print(f"  no .a3m/.csv alignment under {root}. Tree ({len(files)} files):")
        for f in files[:40]:
            print(f"    {f.relative_to(root)}  ({f.stat().st_size}B)")
        if len(files) > 40:
            print(f"    ... and {len(files) - 40} more")
    return None


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
        print("--- boltz stdout (tail) ---")
        print(proc.stdout[-2000:])
        print("--- boltz stderr (tail) ---")
        print(proc.stderr[-3000:])
        print(f"  ! keeping {work} for inspection (not deleted)")
        raise SystemExit(f"boltz predict failed with exit code {proc.returncode}")
    print("  seeding prediction finished; locating the alignment ...", flush=True)

    msa = find_msa(work)
    if msa is None:
        print(f"\n  ! keeping {work} for inspection (not deleted)")
        raise SystemExit(
            f"Could not find a generated MSA under {work}.\n"
            "The tree listing above shows what boltz actually wrote; point\n"
            "Boltz2Config(precomputed_msa=...) at the right file, or re-run with\n"
            "--cache to place it manually. Falling back to --use_msa_server still\n"
            "works but re-queries the server for EVERY prediction.")

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
