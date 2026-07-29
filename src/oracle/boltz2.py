"""
Boltz-2 binding oracle (GPU) — the accuracy edge over AiDTA.

Co-folds a candidate DNA aptamer together with the MMP9 receptor and scores the
predicted INTERFACE (interchain ipTM / chain-pair ipTM), which is a real
protein-nucleic-acid binding signal — replacing AiDTA's secondary-structure
proxy reward.

Boltz-2 is MIT-licensed (commercial-OK) and co-folds protein/DNA/RNA. NOTE its
affinity head is small-molecule-only, so we use the STRUCTURE/interface
confidence, not the affinity module, for nucleic-acid targets. Treat the score
as a strong-but-noisy ranker (per ACS Synth. Biol. 2025) — Stage 4 consensus
(AF3) and Stage 5 MM/GBSA refine it further.

Requires a GPU workstation with Boltz installed:
    pip install boltz            # https://github.com/jwohlwend/boltz
    # first run downloads weights (~cached under ~/.boltz)

This module is intentionally a thin, correct wrapper: it writes Boltz YAML input,
invokes `boltz predict`, and parses the confidence JSON. It raises a clear error
if boltz is unavailable, so it never silently fabricates a score.
"""
from __future__ import annotations
from dataclasses import dataclass
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from oracle.interface import Oracle, OracleScore  # noqa: E402


@dataclass
class Boltz2Config:
    receptor_sequence: str          # MMP9 protein sequence (e.g. catalytic domain)
    protein_id: str = "A"
    dna_id: str = "B"
    ligand_type: str = "dna"     # "dna" or "rna" — the entity Boltz co-folds
    use_msa_server: bool = True     # let boltz build the protein MSA (needs network)
    precomputed_msa: str | None = None   # path to a protein .a3m to reuse (offline)
    devices: int = 1
    out_root: str | None = None     # where boltz writes; temp dir if None
    boltz_bin: str = "boltz"
    extra_args: tuple = ()
    keep_outputs: bool = False   # keep auto-created temp dirs (debugging)


REPO_ROOT = Path(__file__).resolve().parents[2]
MSA_CACHE_DIR = REPO_ROOT / "data" / "mmp9"


def msa_cache_base(receptor_sequence: str) -> Path:
    """Cache path keyed to the receptor SEQUENCE.

    A single fixed filename would hand a 337-aa alignment to a 162-aa receptor
    when the domain is switched, with nothing detecting the mismatch. Resolving
    from the repo root (not the CWD) also keeps the cache usable when the
    pipeline is invoked from elsewhere — otherwise every prediction silently
    re-queries the public MSA server.
    """
    h = hashlib.sha1(receptor_sequence.encode()).hexdigest()[:10]
    return MSA_CACHE_DIR / f"msa_{h}"


def find_cached_msa(receptor_sequence: str):
    """Return a cached alignment for THIS receptor, if one exists."""
    base = msa_cache_base(receptor_sequence)
    for ext in (".a3m", ".csv"):
        p = Path(f"{base}{ext}")
        if p.exists() and p.stat().st_size > 0:
            return str(p.resolve())
    return None


def _write_yaml(path: Path, cfg: Boltz2Config, aptamer_seq: str):
    prot = [f"  - protein:",
            f"      id: {cfg.protein_id}",
            f"      sequence: {cfg.receptor_sequence}"]
    if cfg.precomputed_msa:
        prot.append(f"      msa: {cfg.precomputed_msa}")
    lig = cfg.ligand_type.lower()
    if lig not in ("dna", "rna"):
        raise ValueError(f"ligand_type must be 'dna' or 'rna', got {lig!r}")
    seq = aptamer_seq.replace("&", "").upper()
    # keep the alphabet consistent with the declared entity
    seq = seq.replace("U", "T") if lig == "dna" else seq.replace("T", "U")
    lines = ["version: 1", "sequences:", *prot,
             f"  - {lig}:",
             f"      id: {cfg.dna_id}",
             f"      sequence: {seq}"]
    path.write_text("\n".join(lines) + "\n")


def _find_confidence_json(out_dir: Path, stem: str | None = None):
    """Locate the confidence JSON for a specific prediction.

    Matching on `stem` matters: boltz names outputs after the input file, and a
    directory may hold results for several inputs. Returning the wrong file
    would attribute another sequence's score to this candidate.
    """
    pats = ([f"confidence_{stem}_model_0.json", f"confidence_{stem}*.json"]
            if stem else []) + ["confidence_*_model_0.json", "confidence_*.json"]
    for pat in pats:
        hits = sorted(out_dir.rglob(pat))
        if hits:
            return hits[0]
    return None


class Boltz2Oracle(Oracle):
    def __init__(self, config: Boltz2Config):
        self.cfg = config
        if config.precomputed_msa is None:
            cached = find_cached_msa(config.receptor_sequence)
            if cached:
                self.cfg.precomputed_msa = cached
        if shutil.which(config.boltz_bin) is None:
            raise RuntimeError(
                f"'{config.boltz_bin}' not found on PATH. Install with "
                "`pip install boltz` on the GPU workstation. This oracle "
                "never fabricates scores, so it stops here.")

    def _binding_score(self, conf: dict) -> tuple[float, dict]:
        """Derive a [0,1] interface binding score from Boltz confidence JSON.

        Prefers the protein-DNA chain-pair ipTM (the true interface signal) and
        falls back to the global ipTM. Boltz keys `pair_chains_iptm` by chain
        INDEX as a string ("0"/"1"), not by the chain ids we set in the YAML, so
        we look up both spellings and, failing that, take the single off-diagonal
        entry of a 2-chain matrix. Getting this wrong silently downgrades the
        score to the weaker global ipTM.
        """
        def num(x):
            return float(x) if isinstance(x, (int, float)) else None

        iptm = num(conf.get("iptm"))
        if iptm is None:
            iptm = num(conf.get("complex_iptm"))

        pair = conf.get("pair_chains_iptm") or conf.get("chains_pair_iptm") or {}
        interface = None
        if isinstance(pair, dict) and pair:
            a, b = str(self.cfg.protein_id), str(self.cfg.dna_id)
            for k1, k2 in ((a, b), (b, a), ("0", "1"), ("1", "0")):
                sub = pair.get(k1)
                if isinstance(sub, dict) and num(sub.get(k2)) is not None:
                    interface = num(sub.get(k2))
                    break
            if interface is None:
                # any 2-chain matrix: take the off-diagonal value
                offdiag = [num(v) for k1, sub in pair.items()
                           if isinstance(sub, dict)
                           for k2, v in sub.items()
                           if k1 != k2 and num(v) is not None]
                if offdiag:
                    interface = max(offdiag)

        value = interface if interface is not None else iptm
        value = float(value) if value is not None else 0.0
        return max(0.0, min(1.0, value)), {
            "interface_iptm": interface, "complex_iptm": iptm,
            "used": "chain_pair" if interface is not None else "global_iptm",
            "ptm": conf.get("ptm"), "confidence_score": conf.get("confidence_score"),
        }

    def score(self, sequence: str) -> OracleScore:
        cfg = self.cfg
        auto_dir = cfg.out_root is None
        root = Path(cfg.out_root) if cfg.out_root else Path(tempfile.mkdtemp(prefix="boltz_"))
        root.mkdir(parents=True, exist_ok=True)

        # Name the input after the sequence. Boltz records completed inputs by
        # file stem and SKIPS anything already processed in the same --out_dir
        # ("All inputs are already processed."), exiting 0 without predicting.
        # With a fixed name like input.yaml a re-run into a reused directory
        # would silently return the PREVIOUS sequence's score and structure.
        stem = "apt_" + hashlib.sha1(sequence.replace("&", "").encode()).hexdigest()[:12]
        yaml_path = root / f"{stem}.yaml"
        _write_yaml(yaml_path, cfg, sequence)

        cmd = [cfg.boltz_bin, "predict", str(yaml_path),
               "--out_dir", str(root), "--devices", str(cfg.devices),
               "--output_format", "pdb", "--override"]
        if cfg.use_msa_server and not cfg.precomputed_msa:
            cmd.append("--use_msa_server")
        cmd.extend(cfg.extra_args)

        try:
            proc = subprocess.run(cmd, capture_output=True, text=True)
            if proc.returncode != 0:
                raise RuntimeError(f"boltz predict failed:\n{proc.stderr[-2000:]}")

            conf_json = _find_confidence_json(root, stem)
            if conf_json is None:
                raise RuntimeError(f"No Boltz confidence JSON for {stem} under {root}")
            conf = json.loads(conf_json.read_text())
            value, details = self._binding_score(conf)
            details["confidence_json"] = str(conf_json)
            details["stem"] = stem
            struct = sorted(root.rglob(f"*{stem}*.pdb")) or sorted(root.rglob("*.pdb"))
            details["structure"] = str(struct[0]) if struct else None
            return OracleScore(value=round(value, 4), kind="boltz2_interface",
                               details=details)
        finally:
            # A per-call temp tree holds the MSA, .npz tensors and Lightning logs;
            # a full pilot makes hundreds of them and would fill /tmp on a
            # compute node. Only clean up dirs we created ourselves.
            if auto_dir and cfg.keep_outputs is False:
                shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    # Smoke test of I/O plumbing WITHOUT running boltz (no GPU here).
    cfg = Boltz2Config(receptor_sequence="FVLTEGNPRWEQTHLTYRIENYTPDL")  # placeholder
    p = Path(tempfile.mkdtemp()) / "input.yaml"
    _write_yaml(p, cfg, "TCTGTCTGGGCGACATTTGCCGTGGCCACCCAGCG")
    print("Wrote example Boltz YAML:\n")
    print(p.read_text())
    print("Run on the GPU workstation: boltz predict input.yaml --use_msa_server")
