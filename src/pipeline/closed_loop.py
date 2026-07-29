"""
Closed-loop aptamer design — the core improvement over AiDTA.

AiDTA: generate 20,000 sequences whose reward is fold self-consistency (no
protein), then filter post-hoc. The generator never learns what binds.

Here: generate -> score with an oracle that CONTAINS THE PROTEIN -> keep elites
-> bias the next round toward what the oracle liked. The oracle signal feeds back
into generation, so the search concentrates on binders.

Loop (round r):
  1. Sample candidates by fragment assembly (src/generator/assembler.py) from the
     current fragment weights.
  2. Cheap pre-filter with the CPU structural prior (optional, --prefilter) so the
     expensive oracle only sees plausibly-folded candidates.
  3. Score survivors with the oracle (Boltz-2 co-fold interface score on GPU, or
     the CPU proxy for dry runs).
  4. Keep the top `elite_frac`; update fragment weights by how often each fragment
     appears in elites (cross-entropy / estimation-of-distribution update with
     smoothing) -> next round samples those fragments more often.
  5. Track the global best; stop after --rounds or when improvement stalls.

This is a deliberately simple, robust learner (an EDA / cross-entropy method).
It needs no gradient through the oracle and works with a few hundred oracle calls
per round — the right regime when each call is a GPU co-fold. An MCTS/GFlowNet
generator can be swapped in later behind the same interface.

Usage (GPU server):
  python3 src/pipeline/closed_loop.py \
      --receptor data/mmp9/mmp9_receptor.fasta \
      --pool data/mmp9/pool.csv --rounds 8 --per-round 96 --out results/run1

Dry run on CPU (no binding signal, verifies the machinery):
  python3 src/pipeline/closed_loop.py --oracle proxy --pool data/mmp9/pool.csv \
      --rounds 3 --per-round 40 --out results/dryrun
"""
from __future__ import annotations
import argparse
import json
import random
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from generator.assembler import available, _realized_len  # noqa: E402
from generator.pool import load_pool, build_generic_pool  # noqa: E402
from oracle.interface import StructureProxyOracle  # noqa: E402
from scoring.secondary_structure import score_candidate  # noqa: E402


# --------------------------------------------------------------------------
# weighted assembly
# --------------------------------------------------------------------------
def weighted_assemble(pool, weights, min_len, max_len, rng):
    """Assemble one candidate, choosing fragments in proportion to `weights`.

    Position (linkage site) is chosen uniformly among legal sites, so the learned
    signal lives in the fragment distribution — the part that carries chemistry.
    """
    seq, struct = "", ""
    guard = 0
    while _realized_len(seq) < min_len:
        guard += 1
        if guard > 500:
            return None
        # pick a fragment by weight, retry a few times if it would overshoot
        for _ in range(20):
            idx = rng.choices(range(len(pool)), weights=weights, k=1)[0]
            frag = pool[idx]
            states = available(seq, struct, [frag])
            states = [s for s in states if _realized_len(s[0]) <= max_len]
            if states:
                seq, struct = rng.choice(states)
                break
        else:
            return None
    return seq, struct


def fragments_used(sequence: str, pool) -> set[int]:
    """Which pool fragments appear in the assembled sequence (substring match on
    the realized strands). Approximate but sufficient for a distribution update."""
    flat = sequence.replace("&", "")
    used = set()
    for i, (fs, _) in enumerate(pool):
        core = fs.split("&")[0]
        if len(core) > 1 and core in flat:
            used.add(i)
    return used


# --------------------------------------------------------------------------
# oracle construction
# --------------------------------------------------------------------------
def build_oracle(kind: str, receptor: str | None):
    if kind == "proxy":
        return StructureProxyOracle(), "structure_proxy"
    from oracle.boltz2 import Boltz2Oracle, Boltz2Config
    if not receptor:
        raise SystemExit("--receptor is required for the boltz2 oracle")
    seq = "".join(l.strip() for l in Path(receptor).read_text().splitlines()
                  if l and not l.startswith(">"))
    return Boltz2Oracle(Boltz2Config(receptor_sequence=seq)), "boltz2_interface"


# --------------------------------------------------------------------------
# main loop
# --------------------------------------------------------------------------
def run(args):
    rng = random.Random(args.seed)
    pool = load_pool(args.pool) if args.pool and Path(args.pool).exists() \
        else build_generic_pool()
    oracle, oracle_kind = build_oracle(args.oracle, args.receptor)
    prefilter = StructureProxyOracle() if args.prefilter else None

    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)

    # ---- resume from checkpoint (Slurm allocations are time-boxed) ---------
    ckpt = outdir / "checkpoint.json"
    weights = [1.0] * len(pool)
    seen: dict[str, dict] = {}
    history = []
    start_round = 1
    if args.resume and ckpt.exists():
        state = json.loads(ckpt.read_text())
        if len(state.get("weights", [])) == len(pool):
            weights = state["weights"]
            seen = state.get("seen", {})
            history = state.get("history", [])
            start_round = state.get("next_round", 1)
            print(f"resumed from {ckpt}: round {start_round}, "
                  f"{len(seen)} sequences already scored\n")
        else:
            print(f"!! checkpoint pool size mismatch — starting fresh\n")

    t0 = time.time()

    def save_ckpt(next_round):
        ckpt.write_text(json.dumps(
            {"weights": weights, "seen": seen, "history": history,
             "next_round": next_round, "pool_size": len(pool)},
            default=str))

    print(f"pool={len(pool)} fragments | oracle={oracle_kind} | "
          f"rounds={args.rounds} x {args.per_round}\n")

    for r in range(start_round, args.rounds + 1):
        # ---- 1. sample ----------------------------------------------------
        cands, tries = [], 0
        while len(cands) < args.per_round and tries < args.per_round * 40:
            tries += 1
            res = weighted_assemble(pool, weights, args.min_len, args.max_len, rng)
            if not res:
                continue
            flat = res[0].replace("&", "")
            if flat in seen:
                continue
            cands.append(res)
            seen[flat] = {}
        if not cands:
            print(f"round {r}: no new candidates; stopping")
            break

        # ---- 2. cheap pre-filter -----------------------------------------
        if prefilter and args.prefilter_keep < 1.0:
            scored = [(prefilter.score(s.replace("&", "")).value, (s, st))
                      for s, st in cands]
            scored.sort(key=lambda x: x[0], reverse=True)
            k = max(1, int(len(scored) * args.prefilter_keep))
            cands = [c for _, c in scored[:k]]

        # ---- 3. oracle ----------------------------------------------------
        rows = []
        for s, st in cands:
            flat = s.replace("&", "")
            try:
                sc = oracle.score(flat)
            except Exception as e:
                print(f"  ! oracle failed on {flat[:20]}...: {e}")
                continue
            rec = {"sequence": flat, "assembled": s, "structure": st,
                   "oracle": sc.value, "oracle_kind": sc.kind,
                   "details": sc.details, "round": r}
            rows.append(rec)
            seen[flat] = rec
        if not rows:
            print(f"round {r}: no oracle scores; stopping")
            break

        # ---- 4. elites -> weight update ----------------------------------
        rows.sort(key=lambda x: x["oracle"], reverse=True)
        n_elite = max(1, int(len(rows) * args.elite_frac))
        elites = rows[:n_elite]

        counts = [0.0] * len(pool)
        for e in elites:
            for i in fragments_used(e["assembled"], pool):
                counts[i] += 1.0
        total = sum(counts) or 1.0
        target = [c / total for c in counts]
        uniform = 1.0 / len(pool)
        # smooth toward uniform so no fragment dies out (keeps exploration alive)
        weights = [(1 - args.lr) * w + args.lr * (args.smooth * uniform +
                                                  (1 - args.smooth) * t * len(pool))
                   for w, t in zip(weights, target)]
        weights = [max(w, args.min_weight) for w in weights]

        best = max(seen.values(), key=lambda v: v.get("oracle", -1))
        mean_r = sum(x["oracle"] for x in rows) / len(rows)
        history.append({"round": r, "n_scored": len(rows),
                        "mean": round(mean_r, 4),
                        "elite_mean": round(sum(e["oracle"] for e in elites) / n_elite, 4),
                        "best_so_far": round(best["oracle"], 4)})
        print(f"round {r:2d}: scored={len(rows):3d}  mean={mean_r:.4f}  "
              f"elite={history[-1]['elite_mean']:.4f}  best={best['oracle']:.4f}")
        save_ckpt(r + 1)

        if args.max_hours and (time.time() - t0) / 3600 >= args.max_hours:
            print(f"\nreached --max-hours {args.max_hours}; checkpointed at round {r}. "
                  f"Resume with the same --out and --resume.")
            break

    # ---- report -----------------------------------------------------------
    ranked = sorted((v for v in seen.values() if v.get("oracle") is not None),
                    key=lambda v: v["oracle"], reverse=True)
    for rec in ranked[:args.top]:
        rec["local"] = score_candidate(rec["sequence"])

    result = {
        "oracle": oracle_kind,
        "pool_size": len(pool),
        "rounds": len(history),
        "total_scored": len(ranked),
        "elapsed_sec": round(time.time() - t0, 1),
        "history": history,
        "top": ranked[:args.top],
    }
    (outdir / "result.json").write_text(json.dumps(result, indent=2, default=str))

    with open(outdir / "top_candidates.csv", "w") as f:
        f.write("rank,sequence,length,oracle_score,mfe,gc_fraction\n")
        for i, rec in enumerate(ranked[:args.top], 1):
            loc = rec.get("local", {})
            f.write(f"{i},{rec['sequence']},{len(rec['sequence'])},"
                    f"{rec['oracle']:.4f},{loc.get('mfe','')},"
                    f"{loc.get('gc_fraction','')}\n")

    print(f"\nTop {min(args.top, len(ranked))} candidates -> {outdir}/top_candidates.csv")
    for i, rec in enumerate(ranked[:min(5, len(ranked))], 1):
        print(f"  {i}. {rec['sequence']}  score={rec['oracle']:.4f}")
    if oracle_kind == "structure_proxy":
        print("\nNOTE: dry run with the CPU proxy — these scores contain NO "
              "binding information. Re-run with --oracle boltz2 on the GPU.")
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--oracle", choices=["boltz2", "proxy"], default="boltz2")
    ap.add_argument("--receptor", help="receptor FASTA (required for boltz2)")
    ap.add_argument("--pool", default="data/mmp9/pool.csv")
    ap.add_argument("--rounds", type=int, default=8)
    ap.add_argument("--per-round", type=int, default=96)
    ap.add_argument("--elite-frac", type=float, default=0.2)
    ap.add_argument("--lr", type=float, default=0.5, help="weight update rate")
    ap.add_argument("--smooth", type=float, default=0.3, help="pull toward uniform")
    ap.add_argument("--min-weight", type=float, default=0.05)
    ap.add_argument("--min-len", type=int, default=30)
    ap.add_argument("--max-len", type=int, default=45)
    ap.add_argument("--prefilter", action="store_true",
                    help="cheap CPU structural pre-filter before the oracle")
    ap.add_argument("--prefilter-keep", type=float, default=0.5)
    ap.add_argument("--top", type=int, default=30)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="results/run1")
    ap.add_argument("--resume", action="store_true",
                    help="continue from checkpoint.json in --out (for time-boxed "
                         "Slurm allocations)")
    ap.add_argument("--max-hours", type=float, default=0,
                    help="stop and checkpoint after this many hours (0 = no limit); "
                         "set below your Slurm --time so the run ends cleanly")
    run(ap.parse_args())


if __name__ == "__main__":
    main()
