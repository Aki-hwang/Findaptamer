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
    use_msa_server: bool = True     # let boltz build the protein MSA (needs network)
    precomputed_msa: str | None = None   # path to a protein .a3m to reuse (offline)
    devices: int = 1
    out_root: str | None = None     # where boltz writes; temp dir if None
    boltz_bin: str = "boltz"
    extra_args: tuple = ()


DEFAULT_MSA_CACHE = "data/mmp9/mmp9_receptor_msa"


def find_cached_msa(base: str = DEFAULT_MSA_CACHE):
    """Return a cached receptor alignment (.a3m/.csv) if one exists.

    Generated once by src/target/make_msa.py. Reusing it removes the per-call
    MSA-server round trip, which otherwise dominates runtime in the closed loop.
    """
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
    lines = ["version: 1", "sequences:", *prot,
             f"  - dna:",
             f"      id: {cfg.dna_id}",
             f"      sequence: {aptamer_seq.replace('&', '')}"]
    path.write_text("\n".join(lines) + "\n")


def _find_confidence_json(out_dir: Path):
    hits = list(out_dir.rglob("confidence_*_model_0.json"))
    if not hits:
        hits = list(out_dir.rglob("confidence_*.json"))
    return hits[0] if hits else None


class Boltz2Oracle(Oracle):
    def __init__(self, config: Boltz2Config):
        self.cfg = config
        if config.precomputed_msa is None:
            cached = find_cached_msa()
            if cached:
                self.cfg.precomputed_msa = cached
        if shutil.which(config.boltz_bin) is None:
            raise RuntimeError(
                f"'{config.boltz_bin}' not found on PATH. Install with "
                "`pip install boltz` on the GPU workstation. This oracle "
                "never fabricates scores, so it stops here.")

    def _binding_score(self, conf: dict) -> tuple[float, dict]:
        """Derive a [0,1] interface binding score from Boltz confidence JSON.

        Prefer the protein-DNA chain-pair ipTM (true interface signal); fall back
        to global ipTM. Also surface complex ipTM/ptm for transparency.
        """
        iptm = conf.get("iptm") or conf.get("complex_iptm")
        pair = conf.get("pair_chains_iptm") or conf.get("chains_pair_iptm") or {}
        # pair is typically a nested dict {chainA: {chainB: value}}
        interface = None
        try:
            a, b = self.cfg.protein_id, self.cfg.dna_id
            interface = (pair.get(a, {}).get(b)
                         if isinstance(pair, dict) else None)
            if interface is None and isinstance(pair, dict):
                interface = pair.get(b, {}).get(a)
        except Exception:
            interface = None
        value = interface if interface is not None else iptm
        value = float(value) if value is not None else 0.0
        return max(0.0, min(1.0, value)), {
            "interface_iptm": interface, "complex_iptm": iptm,
            "ptm": conf.get("ptm"), "confidence_score": conf.get("confidence_score"),
        }

    def score(self, sequence: str) -> OracleScore:
        cfg = self.cfg
        root = Path(cfg.out_root) if cfg.out_root else Path(tempfile.mkdtemp(prefix="boltz_"))
        root.mkdir(parents=True, exist_ok=True)
        yaml_path = root / "input.yaml"
        _write_yaml(yaml_path, cfg, sequence)

        cmd = [cfg.boltz_bin, "predict", str(yaml_path),
               "--out_dir", str(root), "--devices", str(cfg.devices),
               "--output_format", "pdb"]
        if cfg.use_msa_server and not cfg.precomputed_msa:
            cmd.append("--use_msa_server")
        cmd.extend(cfg.extra_args)

        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"boltz predict failed:\n{proc.stderr[-2000:]}")

        conf_json = _find_confidence_json(root)
        if conf_json is None:
            raise RuntimeError(f"No Boltz confidence JSON under {root}")
        conf = json.loads(conf_json.read_text())
        value, details = self._binding_score(conf)
        details["confidence_json"] = str(conf_json)
        return OracleScore(value=round(value, 4), kind="boltz2_interface", details=details)


if __name__ == "__main__":
    # Smoke test of I/O plumbing WITHOUT running boltz (no GPU here).
    cfg = Boltz2Config(receptor_sequence="FVLTEGNPRWEQTHLTYRIENYTPDL")  # placeholder
    p = Path(tempfile.mkdtemp()) / "input.yaml"
    _write_yaml(p, cfg, "TCTGTCTGGGCGACATTTGCCGTGGCCACCCAGCG")
    print("Wrote example Boltz YAML:\n")
    print(p.read_text())
    print("Run on the GPU workstation: boltz predict input.yaml --use_msa_server")
