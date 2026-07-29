"""
Fragment-assembly engine for aptamer design.

Faithful, self-contained port of the AiDTA assembly semantics (game.py
`available`) with the reward re-implemented on ViennaRNA (DNA parameters) via
`src/scoring/secondary_structure.py`. No RNAstructure or GPU dependency, so it
runs on CPU.

State representation (same as AiDTA):
  A partially assembled aptamer is a pair (sequence, structure) of '&'-joined
  strand segments, e.g.
      sequence  = "GGC&GCC"          structure = "(((&)))"
  The '&' marks strand breaks introduced by double-stranded fragments (the two
  arms of a stem). The realized DNA sequence is `sequence.replace('&','')`.

Actions:
  Insert a fragment (its sequence+structure) at the end of one of the current
  segments. This reproduces AiDTA's `available()` exactly, including how a ds
  fragment adds a new '&' break.

Provides:
  * available(sequence, structure, pool)      -> list of (seq, struct) next states
  * random_assemble(pool, min_len, ...)        -> AiDTA's random-assembly control
  * greedy_assemble(pool, min_len, ...)        -> structure-consistency-greedy
  * batch_generate(pool, n, strategy, ...)     -> many candidates + local scores

The closed-loop driver (src/pipeline/closed_loop.py) reuses `available()` as its
action space and replaces the greedy structural criterion with an oracle that
contains the protein.
"""
from __future__ import annotations
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scoring.secondary_structure import fold_dna, structure_similarity  # noqa: E402


def available(sequence: str, structure: str, pool):
    """Enumerate legal next states by inserting each fragment at each segment end.

    `pool` is a list of (fragment_sequence, fragment_structure). Returns a list of
    (new_sequence, new_structure). Faithful to AiDTA game.py `available()`.
    """
    results = []
    seq_parts = sequence.split("&")
    struc_parts = structure.split("&")
    for frag_seq, frag_struct in pool:
        for idx in range(len(seq_parts)):
            new_seq_parts = seq_parts.copy()
            new_struc_parts = struc_parts.copy()
            new_seq_parts[idx] = new_seq_parts[idx] + frag_seq
            new_struc_parts[idx] = new_struc_parts[idx] + frag_struct
            results.append(("&".join(new_seq_parts), "&".join(new_struc_parts)))
    return results


def _realized_len(sequence: str) -> int:
    return len(sequence.replace("&", ""))


def random_assemble(pool, min_len=50, max_len=61, rng=None):
    """AiDTA's random-assembly control: pick random legal insertions until the
    realized length reaches [min_len, max_len]. Returns (sequence, structure) or
    None if it cannot finish within the length budget."""
    rng = rng or random
    seq, struct = "", ""
    guard = 0
    while _realized_len(seq) < min_len:
        guard += 1
        if guard > 2000:
            return None
        states = available(seq, struct, pool)
        states = [s for s in states if _realized_len(s[0]) <= max_len]
        if not states:
            return None
        seq, struct = rng.choice(states)
    return seq, struct


def greedy_assemble(pool, min_len=50, max_len=61, beam=8, rng=None):
    """Structure-consistency-greedy assembly.

    At each step, among candidate insertions, keep the one whose folded DNA
    structure best agrees (paired/unpaired) with the assembled fragment structure
    — i.e. AiDTA's reward, used greedily instead of via RL. `beam` samples a
    random subset of actions per step to keep it cheap and diverse.

    NOTE: this criterion is the one we show to be degenerate (an all-unpaired
    sequence scores 1.0). It is kept as a baseline/control, not as the objective.
    """
    rng = rng or random
    seq, struct = "", ""
    guard = 0
    while _realized_len(seq) < min_len:
        guard += 1
        if guard > 500:
            return None
        states = [s for s in available(seq, struct, pool)
                  if _realized_len(s[0]) <= max_len]
        if not states:
            return None
        if len(states) > beam:
            states = rng.sample(states, beam)
        best, best_sim = None, -1.0
        for cand_seq, cand_struct in states:
            folded, _ = fold_dna(cand_seq)
            sim = structure_similarity(folded, cand_struct)
            if sim > best_sim:
                best, best_sim = (cand_seq, cand_struct), sim
        seq, struct = best
    return seq, struct


def batch_generate(pool, n=100, strategy="greedy", min_len=50, max_len=61,
                   seed=0, dedup=True):
    """Generate n candidates and attach local scores. Returns list of dicts."""
    from scoring.secondary_structure import score_candidate
    rng = random.Random(seed)
    out, seen = [], set()
    attempts = 0
    while len(out) < n and attempts < n * 50:
        attempts += 1
        res = (greedy_assemble if strategy == "greedy" else random_assemble)(
            pool, min_len, max_len, rng=rng)
        if res is None:
            continue
        seq, struct = res
        flat = seq.replace("&", "")
        if dedup and flat in seen:
            continue
        seen.add(flat)
        sc = score_candidate(seq, assembled_structure=struct)
        sc["sequence"] = flat
        sc["assembled_structure"] = struct.replace("&", "")
        sc["strategy"] = strategy
        out.append(sc)
    return out


if __name__ == "__main__":
    pool = [
        ("A", "."), ("G", "."), ("C", "."), ("T", "."),
        ("GGC&GCC", "(((&)))"), ("CAC&GTG", "(((&)))"),
        ("GTGG&CCAC", "((((&))))"), ("GGGTG&CACCC", "(((((&)))))"),
        ("TCTG", "...."), ("GAGA", "...."), ("TTTTAC", "......"),
    ]
    cands = batch_generate(pool, n=5, strategy="greedy", seed=1)
    for c in cands:
        print(f"{c['sequence']}  sim={c['assembly_similarity']:.2f}  "
              f"mfe={c['mfe']:.1f}  ensfreq={c['ens_freq_mfe']:.2f}")
