"""
Chai-1 co-folding oracle — the independent second model for Stage-4 consensus.

Why a second model: co-fold confidence is a strong-but-noisy aptamer screener
(ACS Synth. Biol. 2025 found AF3 models aptamer structures well yet "fails to
identify promising aptamers in screening cases"). Agreement between two
*independently trained* predictors is the cheapest available antidote to
single-model overfitting, and it is what Stage 4 tests.

Chai-1 is Apache-2.0 (weights included) and predicts protein/DNA/RNA complexes
without nucleic-acid MSAs. Install on the GPU node:

    pip install chai_lab

Like Boltz2Oracle, this never fabricates a score: if chai_lab is missing or the
run fails, it raises.
"""
from __future__ import annotations
from dataclasses import dataclass
import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from oracle.interface import Oracle, OracleScore  # noqa: E402


@dataclass
class Chai1Config:
    receptor_sequence: str
    out_root: str | None = None
    num_trunk_recycles: int = 3
    num_diffn_timesteps: int = 200
    seed: int = 0
    device: str = "cuda:0"
    use_esm_embeddings: bool = True


def write_chai_fasta(path: Path, receptor: str, aptamer: str):
    """Chai-1 takes a FASTA whose headers declare the entity type."""
    path.write_text(
        f">protein|name=receptor\n{receptor}\n"
        f">dna|name=aptamer\n{aptamer.replace('&', '')}\n"
    )


class Chai1Oracle(Oracle):
    def __init__(self, config: Chai1Config):
        self.cfg = config
        try:
            import chai_lab  # noqa: F401
        except ImportError as e:
            raise RuntimeError(
                "chai_lab is not installed. On the GPU node run "
                "`pip install chai_lab`. This oracle never fabricates scores, "
                f"so it stops here. ({e})")

    def score(self, sequence: str) -> OracleScore:
        from chai_lab.chai1 import run_inference
        import torch

        cfg = self.cfg
        root = Path(cfg.out_root) if cfg.out_root else Path(
            tempfile.mkdtemp(prefix="chai_"))
        root.mkdir(parents=True, exist_ok=True)
        fasta = root / "input.fasta"
        write_chai_fasta(fasta, cfg.receptor_sequence, sequence)

        candidates = run_inference(
            fasta_file=fasta,
            output_dir=root,
            num_trunk_recycles=cfg.num_trunk_recycles,
            num_diffn_timesteps=cfg.num_diffn_timesteps,
            seed=cfg.seed,
            device=torch.device(cfg.device),
            use_esm_embeddings=cfg.use_esm_embeddings,
        )

        # ranking data carries aggregate + per-chain-pair scores
        best_iptm, details = 0.0, {}
        try:
            scores = candidates.ranking_data[0]
            agg = getattr(scores, "aggregate_score", None)
            iptm_obj = getattr(scores, "ptm_scores", None)
            best_iptm = float(getattr(iptm_obj, "interface_ptm", agg))
            details = {"aggregate_score": float(agg) if agg is not None else None,
                       "interface_ptm": best_iptm}
        except Exception:
            # fall back to the scores json Chai writes alongside the structures
            for js in sorted(root.rglob("scores*.json")):
                d = json.loads(js.read_text())
                best_iptm = float(d.get("iptm") or d.get("interface_ptm") or 0.0)
                details = d
                break

        cif = sorted(root.rglob("*.cif"))
        if cif:
            details["structure"] = str(cif[0])
        return OracleScore(value=max(0.0, min(1.0, best_iptm)),
                           kind="chai1_interface", details=details)


if __name__ == "__main__":
    p = Path(tempfile.mkdtemp()) / "input.fasta"
    write_chai_fasta(p, "FVLTEGNPRWEQTHLTYRIENYTPDL", "TCTGTCTGGGCGACATTTGCC")
    print("Example Chai-1 FASTA:\n")
    print(p.read_text())
    print("Run on the GPU node after: pip install chai_lab")
