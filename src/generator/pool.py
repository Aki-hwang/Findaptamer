"""
Fragment-pool construction for the assembler (Stage 1 front end).

The 10,880-fragment library (`src/fragments/build_library.py`) is
target-independent. AiDTA turns it into a small *target pool* (~44 fragments in
the paper) by docking every fragment to the protein epitope with HDOCK and
keeping the top scorers; the generator then assembles aptamers out of that pool.

That docking step needs external servers/binaries and is the one place the paper
injects target information *before* generation. We support two ways to get a pool
so the closed loop can run with or without a GPU/docking service:

  * build_generic_pool()  — CPU, no docking. Selects a compact, structurally
    diverse set of ss + ds fragments straight from the library, ranked by a
    cheap intrinsic-foldability score. This lets the ORACLE do the targeting in
    the loop (our thesis: the in-loop binding oracle, not pre-docking, is what
    should steer sequence toward the epitope). Deterministic (seeded).

  * load_docked_pool(csv) — parse an HDOCK/HADDOCK results table
    (fragment, docking_score[, ...]) and keep the top-N by score. Use this on the
    workstation once fragments have been docked to the MMP9 epitope, to reproduce
    the AiDTA-style target-biased pool. No scores are ever fabricated: if the
    file is missing/empty this raises rather than returning a guessed pool.

A Fragment is a small immutable record; a Pool is just a list of Fragments plus
provenance so downstream code can report how the pool was built.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from fragments.build_library import build  # noqa: E402
from scoring.secondary_structure import fold_dna  # noqa: E402


@dataclass(frozen=True)
class Fragment:
    sequence: str          # DNA, '&'-joined for duplexes (e.g. "GGC&GCC")
    structure: str         # dot-bracket in the same '&'-joined notation
    kind: str              # "ss" or "ds"
    score: float = 0.0     # provenance score (docking or intrinsic); higher = better
    source: str = ""       # how this fragment entered the pool

    @property
    def strands(self) -> list[str]:
        return self.sequence.split("&")

    @property
    def n_nt(self) -> int:
        return sum(len(s) for s in self.strands)


@dataclass
class Pool:
    fragments: list[Fragment]
    provenance: dict = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.fragments)

    def by_kind(self, kind: str) -> list[Fragment]:
        return [f for f in self.fragments if f.kind == kind]

    def summary(self) -> dict:
        return {
            "n": len(self.fragments),
            "n_ss": len(self.by_kind("ss")),
            "n_ds": len(self.by_kind("ds")),
            "n_nt_range": [min((f.n_nt for f in self.fragments), default=0),
                           max((f.n_nt for f in self.fragments), default=0)],
            **self.provenance,
        }


def _gc(seq: str) -> float:
    s = seq.replace("&", "")
    return (s.count("G") + s.count("C")) / len(s) if s else 0.0


def _intrinsic_ss_score(seq: str) -> float:
    """Cheap foldability prior for a single-stranded fragment.

    We want ss fragments that are NOT self-complementary hairpins (those belong
    in the ds set) yet have balanced composition and no homopolymer runs, so the
    assembler can use them as loops/linkers that actually stay single-stranded.
    Higher = better. Purely intrinsic (no target), deterministic.
    """
    s = seq.replace("&", "")
    if not s:
        return 0.0
    ss, mfe = fold_dna(s)
    unpaired_frac = ss.count(".") / len(ss)
    gc = _gc(s)
    gc_balance = 1.0 - abs(gc - 0.5) * 2.0          # 1 at 50% GC, 0 at 0/100%
    # penalize homopolymer runs (synthesis + folding pathology)
    longest_run = 1
    run = 1
    for a, b in zip(s, s[1:]):
        run = run + 1 if a == b else 1
        longest_run = max(longest_run, run)
    run_penalty = 1.0 if longest_run <= 2 else 1.0 / (longest_run - 1)
    return round(0.6 * unpaired_frac + 0.25 * gc_balance + 0.15 * run_penalty, 4)


def _intrinsic_ds_score(seq: str) -> float:
    """Cheap stability prior for a duplex fragment: reward a clean, stable
    Watson-Crick stem with moderate-to-high GC and no wobble ambiguity.
    """
    s1 = seq.split("&")[0]
    if not s1:
        return 0.0
    gc = _gc(s1)
    gc_stability = min(gc / 0.6, 1.0)               # saturates at 60% GC
    longest_run = 1
    run = 1
    for a, b in zip(s1, s1[1:]):
        run = run + 1 if a == b else 1
        longest_run = max(longest_run, run)
    run_penalty = 1.0 if longest_run <= 2 else 1.0 / (longest_run - 1)
    return round(0.8 * gc_stability + 0.2 * run_penalty, 4)


def build_generic_pool(n_ss: int = 30, n_ds: int = 14,
                       lengths=(4, 5, 6)) -> Pool:
    """Compact, structurally sensible pool selected from the full library (CPU).

    Ranks ss and ds fragments of the requested lengths by an intrinsic
    foldability/stability score and keeps the top n_ss / n_ds. Ties are broken
    lexicographically so the result is deterministic (no RNG). The default
    n_ss+n_ds = 44 mirrors the paper's target-pool size, but here the pool is
    target-AGNOSTIC on purpose: targeting is delegated to the in-loop oracle.
    """
    ss_all, ds_all = build(lengths=lengths)
    ss_scored = [Fragment(seq, st, "ss", _intrinsic_ss_score(seq), "generic:intrinsic")
                 for seq, st in ss_all]
    ds_scored = [Fragment(seq, st, "ds", _intrinsic_ds_score(seq), "generic:intrinsic")
                 for seq, st in ds_all]
    ss_top = sorted(ss_scored, key=lambda f: (-f.score, f.sequence))[:n_ss]
    ds_top = sorted(ds_scored, key=lambda f: (-f.score, f.sequence))[:n_ds]
    frags = ss_top + ds_top
    return Pool(frags, provenance={
        "method": "generic_intrinsic",
        "target_biased": False,
        "lengths": list(lengths),
        "note": "targeting delegated to the in-loop oracle, not to pre-docking",
    })


def load_docked_pool(results_csv: str, top_n: int = 44,
                     seq_col: str = "sequence", struct_col: str = "structure",
                     score_col: str = "docking_score",
                     kind_col: str = "type",
                     higher_is_better: bool = True) -> Pool:
    """Load an AiDTA-style target pool from a docking results table.

    Expects a CSV with at least sequence, structure, type, and a docking-score
    column (HDOCK scores are negative and more-negative = better, so pass
    higher_is_better=False for raw HDOCK). Keeps the top_n by score. Raises if
    the file is missing or has no usable rows — it never invents docking scores.
    """
    path = Path(results_csv)
    if not path.exists():
        raise FileNotFoundError(
            f"Docked-pool file not found: {path}. Dock the fragment library to "
            "the MMP9 epitope (HDOCK/HADDOCK) on the workstation first, or use "
            "build_generic_pool() to run the oracle-targeted loop without docking.")
    rows = list(csv.DictReader(path.open()))
    frags = []
    for r in rows:
        if seq_col not in r or score_col not in r:
            continue
        try:
            raw = float(r[score_col])
        except (TypeError, ValueError):
            continue
        seq = r[seq_col].strip()
        st = r.get(struct_col, "").strip() or _implied_structure(seq)
        kind = r.get(kind_col, "").strip() or ("ds" if "&" in seq else "ss")
        score = raw if higher_is_better else -raw
        frags.append(Fragment(seq, st, kind, round(score, 4), f"docked:{path.name}"))
    if not frags:
        raise ValueError(f"No usable rows in {path} (need '{seq_col}' and "
                         f"'{score_col}' columns).")
    frags.sort(key=lambda f: (-f.score, f.sequence))
    frags = frags[:top_n]
    return Pool(frags, provenance={
        "method": "hdock_docked",
        "target_biased": True,
        "source_file": str(path),
        "top_n": top_n,
    })


def _implied_structure(seq: str) -> str:
    """Structure implied by the '&' notation when a table omits it: a duplex is
    fully paired, a single strand is fully unpaired."""
    if "&" in seq:
        a, b = seq.split("&", 1)
        return f"{'(' * len(a)}&{')' * len(b)}"
    return "." * len(seq)


if __name__ == "__main__":
    import json
    pool = build_generic_pool()
    print(json.dumps(pool.summary(), indent=2))
    print("\nTop 5 ss fragments:")
    for f in pool.by_kind("ss")[:5]:
        print(f"  {f.sequence:8s} {f.structure:8s} score={f.score}")
    print("Top 5 ds fragments:")
    for f in pool.by_kind("ds")[:5]:
        print(f"  {f.sequence:14s} {f.structure:10s} score={f.score}")
