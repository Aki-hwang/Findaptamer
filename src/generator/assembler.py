"""
Faithful AiDTA-style assembly engine (CPU) — and the demonstration of *why* its
reward is not enough.

AiDTA (Guo et al., bioRxiv 2025) builds a single-stranded DNA aptamer by
*assembling* fragments from a target pool with an RL policy (MCTS/AlphaZero),
using a reward that is purely SECONDARY-STRUCTURE SELF-CONSISTENCY: fold the
assembled sequence and compare it to the structure implied by concatenating the
constituent fragments' own structures; if the agreement is >= 0.9 the assembly
is "successful" (reward +1). The protein is never seen during assembly.

This module reproduces that assembly action space faithfully:

  * ss fragment  -> appended as an unpaired stretch  ("...")
  * ds fragment  -> opens an intramolecular STEM: strand1 goes down now as
                    "(((", and its reverse-complement partner strand2 is pushed
                    onto a stack to be laid down later as ")))" (a hairpin). ss
                    fragments placed in between become the loop.

The assembled aptamer is therefore a linear ssDNA whose *implied* structure is a
set of hairpins/loops. `assembly_reward` (from src/scoring) then folds the
sequence for real and scores agreement — exactly AiDTA's reward.

THE POINT: run `python3 src/generator/assembler.py`. An assembler that only ever
appends ss fragments produces unstructured sequences whose implied structure is
all-"." — which the folder also predicts as mostly unpaired — so they collect
AiDTA reward ~+1 while being non-binding, unstructured junk. The same sequences
score LOW under `StructureProxyOracle`, and would score low under the Boltz-2
binding oracle. That gap is the degeneracy our oracle-in-the-loop closes.

The engine is deterministic given a seeded `random.Random`, so the closed-loop
driver can reproduce and evolve populations.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from generator.pool import Pool, Fragment, build_generic_pool  # noqa: E402
from scoring.secondary_structure import assembly_reward, fold_dna  # noqa: E402


@dataclass
class Candidate:
    sequence: str                       # linear ssDNA (no '&')
    implied_structure: str              # dot-bracket implied by the fragments
    trace: list = field(default_factory=list)   # ordered (action, fragment) log

    @property
    def length(self) -> int:
        return len(self.sequence)


@dataclass
class AssemblyConfig:
    length_range: tuple = (25, 45)      # target aptamer length (nt)
    max_stems: int = 4                  # max simultaneously-open hairpin stems
    ss_only: bool = False               # if True, never open stems (degeneracy demo)
    # action propensities (renormalized against what is legal at each step)
    p_open: float = 0.40
    p_close: float = 0.30
    p_ss: float = 0.30


class Assembler:
    """Assembles aptamers from a Pool. Oracle-agnostic: it only produces
    candidates + their AiDTA implied structure; scoring is someone else's job."""

    def __init__(self, pool: Pool | None = None, config: AssemblyConfig | None = None):
        self.pool = pool or build_generic_pool()
        self.cfg = config or AssemblyConfig()
        self._ss = self.pool.by_kind("ss")
        self._ds = self.pool.by_kind("ds")
        if not self._ss:
            raise ValueError("Pool has no single-stranded fragments to build loops.")
        if not self._ds and not self.cfg.ss_only:
            # structured assembly needs stems; fall back to ss-only but flag it
            self.cfg.ss_only = True

    # -- assembly ------------------------------------------------------------
    def assemble(self, rng: random.Random) -> Candidate:
        cfg = self.cfg
        seq: list[str] = []
        struct: list[str] = []
        stack: list[tuple[str, str]] = []   # (closing_strand, closing_brackets)
        trace: list = []
        target_len = rng.randint(*cfg.length_range)

        def cur_len() -> int:
            return sum(len(s) for s in seq)

        # reserve room to close every open stem before hitting target_len
        def reserved() -> int:
            return sum(len(cs) for cs, _ in stack)

        while True:
            n = cur_len()
            # must stop once we've reached the target (then close remaining stems)
            if n >= target_len:
                break
            room = target_len - n - reserved()
            legal = []
            weights = []
            # append ss (loop/linker) — legal if a fragment fits in remaining room
            ss_fit = [f for f in self._ss if len(f.sequence) <= max(room, 1)]
            if ss_fit:
                legal.append("ss"); weights.append(cfg.p_ss)
            # open a new stem — legal if allowed, not ss_only, room for stem+partner
            if (not cfg.ss_only and len(stack) < cfg.max_stems and self._ds):
                ds_fit = [f for f in self._ds
                          if 2 * len(f.strands[0]) <= max(room, 0)]
                if ds_fit:
                    legal.append(("open", ds_fit)); weights.append(cfg.p_open)
            # close the innermost stem — legal if one is open and it fits
            if stack and len(stack[-1][0]) <= max(room + reserved(), 1):
                legal.append("close"); weights.append(cfg.p_close)

            if not legal:
                break
            choice = rng.choices(legal, weights=weights, k=1)[0]

            if choice == "ss":
                f = rng.choice(ss_fit)
                seq.append(f.sequence); struct.append("." * len(f.sequence))
                trace.append(("ss", f.sequence))
            elif isinstance(choice, tuple) and choice[0] == "open":
                f = rng.choice(choice[1])
                s1, s2 = f.strands[0], f.strands[1]
                seq.append(s1); struct.append("(" * len(s1))
                stack.append((s2, ")" * len(s2)))
                trace.append(("open", f.sequence))
            elif choice == "close":
                cs, cb = stack.pop()
                seq.append(cs); struct.append(cb)
                trace.append(("close", cs))

        # close any stems still open, so the structure is balanced
        while stack:
            cs, cb = stack.pop()
            seq.append(cs); struct.append(cb)
            trace.append(("close", cs))

        return Candidate("".join(seq), "".join(struct), trace)

    def batch_generate(self, n: int, rng: random.Random) -> list[Candidate]:
        return [self.assemble(rng) for _ in range(n)]


def score_assembly(cand: Candidate, threshold: float = 0.9) -> dict:
    """AiDTA reward for one candidate: fold for real, compare to implied."""
    reward, sim, folded = assembly_reward(cand.sequence, cand.implied_structure,
                                          threshold=threshold)
    return {
        "sequence": cand.sequence,
        "length": cand.length,
        "implied_structure": cand.implied_structure,
        "folded_structure": folded,
        "aidta_similarity": round(sim, 4),
        "aidta_reward": reward,          # +1 "successful assembly" / -1 otherwise
    }


def _demo():
    """Show the degeneracy: ss-only assemblies collect AiDTA reward +1 while
    being unstructured, whereas stem-using assemblies are genuinely folded."""
    import sys as _sys
    _sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from oracle.interface import StructureProxyOracle

    rng = random.Random(7)
    pool = build_generic_pool()
    oracle = StructureProxyOracle()

    print("=== AiDTA reward is gameable: ss-only ('unstructured') assemblies ===")
    ss_only = Assembler(pool, AssemblyConfig(ss_only=True, length_range=(30, 40)))
    passers = []          # oracle values for candidates AiDTA calls "successful"
    worst = None          # airtight case: AiDTA sim ~1.0 but poor structural fitness
    for _ in range(400):
        c = ss_only.assemble(rng)
        a = score_assembly(c)
        ov = oracle.score(c.sequence).value
        if a["aidta_reward"] == 1:               # AiDTA would accept this aptamer
            passers.append(ov)
            if worst is None or a["aidta_similarity"] > worst[1]:
                worst = (a, a["aidta_similarity"], ov)
    if passers:
        print(f"  of the ss-only assemblies AiDTA ACCEPTS (reward +1, n={len(passers)}):")
        print(f"    mean StructureProxyOracle value = {sum(passers)/len(passers):.3f}"
              f"  (min {min(passers):.3f})")
        print("    -> passing AiDTA's structure-consistency reward does NOT imply")
        print("       structural fitness; these are unstructured non-binders.")
    if worst:
        a = worst[0]
        print(f"  airtight example (AiDTA similarity={a['aidta_similarity']}, "
              f"oracle={worst[2]:.3f}):")
        print(f"    seq    : {a['sequence']}")
        print(f"    folded : {a['folded_structure']}  (implied all-unpaired)")

    print("\n=== stem-using assemblies: genuinely structured ===")
    structured = Assembler(pool, AssemblyConfig(ss_only=False, length_range=(30, 40)))
    rewards2, oracle_vals2 = [], []
    for _ in range(200):
        c = structured.assemble(rng)
        a = score_assembly(c)
        rewards2.append(a["aidta_reward"])
        oracle_vals2.append(oracle.score(c.sequence).value)
    succ2 = sum(1 for r in rewards2 if r == 1)
    print(f"  stems : AiDTA 'successful assembly' rate = {succ2}/200 = {succ2/200:.0%}")
    print(f"          mean StructureProxyOracle value  = {sum(oracle_vals2)/len(oracle_vals2):.3f}")

    print("\n  Example stem-using candidate:")
    c = structured.assemble(rng)
    a = score_assembly(c)
    print(f"    seq     : {a['sequence']}")
    print(f"    implied : {a['implied_structure']}")
    print(f"    folded  : {a['folded_structure']}")
    print(f"    aidta_similarity={a['aidta_similarity']}  oracle={oracle.score(c.sequence).value}")


if __name__ == "__main__":
    _demo()
