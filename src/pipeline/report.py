"""
Final report — the synthesis shortlist.

Merges every stage into one ranked table and writes a human-readable report:

  pilot      oracle interface score (Boltz-2, in the generation loop)
  consensus  multi-seed / multi-model agreement + contact-set stability
  energy     MM/GBSA dG ranking

Final ranking is by the physics tier where available (most calibrated), then
consensus, then pilot — i.e. the most trustworthy signal a candidate reached.
Also reports the local DNA-fold properties needed to judge synthesizability and
compares against the literature benchmark (F3B, Kd ~20 nM).

Usage:
  python src/pipeline/report.py --out results/report --n 20
"""
from __future__ import annotations
import argparse
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scoring.secondary_structure import score_candidate  # noqa: E402
from target.mmp9 import BENCHMARK_APTAMERS  # noqa: E402


def load_json(p: str):
    path = Path(p)
    return json.loads(path.read_text()) if path.exists() else None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pilot", default="results/pilot/result.json")
    ap.add_argument("--consensus", default="results/consensus/consensus.json")
    ap.add_argument("--energy", default="results/energy/energy.json")
    ap.add_argument("--n", type=int, default=20, help="shortlist size")
    ap.add_argument("--out", default="results/report")
    args = ap.parse_args()

    pilot = load_json(args.pilot)
    consensus = load_json(args.consensus)
    energy = load_json(args.energy)

    merged: dict[str, dict] = {}
    if pilot:
        for rec in pilot.get("top", []):
            merged.setdefault(rec["sequence"], {})["pilot_score"] = rec.get("oracle")
    if consensus:
        for rec in consensus:
            m = merged.setdefault(rec["sequence"], {})
            m.update({"consensus_score": rec.get("consensus_score"),
                      "score_mean": rec.get("score_mean"),
                      "score_std": rec.get("score_std"),
                      "contact_jaccard": rec.get("contact_jaccard"),
                      "n_contacts": rec.get("n_contacts"),
                      "fnII_fraction": rec.get("fnII_fraction"),
                      "consensus_pass": rec.get("pass")})
    if energy:
        for rec in energy:
            m = merged.setdefault(rec["sequence"], {})
            m.update({"dG_bind_kcal_mol": rec.get("dG_bind_kcal_mol"),
                      "dG_sem": rec.get("sem")})

    if not merged:
        raise SystemExit("no stage outputs found — run the pilot first")

    rows = []
    for seq, d in merged.items():
        loc = score_candidate(seq)
        # tier: 3 = has energy, 2 = has consensus, 1 = pilot only
        tier = 3 if d.get("dG_bind_kcal_mol") is not None else (
            2 if d.get("consensus_score") is not None else 1)
        rows.append({"sequence": seq, "length": len(seq), "tier": tier,
                     **d,
                     "mfe": round(loc["mfe"], 2),
                     "gc_fraction": round(loc["gc_fraction"], 3),
                     "ens_freq_mfe": round(loc["ens_freq_mfe"], 3),
                     "mfe_structure": loc["mfe_structure"]})

    def _num(r, key, default=0.0):
        """None-safe numeric access: a present-but-null value must not reach the
        sort comparison (float < None raises), which would kill the final report
        after every expensive stage had already succeeded."""
        v = r.get(key)
        return default if v is None else float(v)

    def sort_key(r):
        return (-r["tier"],
                _num(r, "dG_bind_kcal_mol"),             # more negative first
                -_num(r, "consensus_score"),
                -_num(r, "pilot_score"))
    rows.sort(key=sort_key)
    short = rows[:args.n]

    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "shortlist.json").write_text(json.dumps(short, indent=2, default=str))

    cols = ["rank", "sequence", "length", "tier", "dG_bind_kcal_mol", "dG_sem",
            "consensus_score", "score_mean", "score_std", "contact_jaccard",
            "fnII_fraction", "pilot_score", "mfe", "gc_fraction", "ens_freq_mfe"]
    with open(outdir / "shortlist.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for i, r in enumerate(short, 1):
            w.writerow([i] + [r.get(c, "") for c in cols[1:]])

    # ---- readable report --------------------------------------------------
    bench = "\n".join(
        f"  - {a['name']:14s} Kd={str(a['affinity_nM'] or '?'):>6s} nM  "
        f"[{a['type']}]  {a['function']}"
        for a in BENCHMARK_APTAMERS)
    stage_note = []
    if pilot:
        stage_note.append(f"pilot: {pilot.get('total_scored')} sequences scored, "
                          f"{pilot.get('rounds')} rounds, oracle={pilot.get('oracle')}")
    if consensus:
        stage_note.append(f"consensus: {sum(1 for r in consensus if r.get('pass'))}"
                          f"/{len(consensus)} passed")
    if energy:
        stage_note.append(f"energy: {len(energy)} candidates with MM/GBSA dG")

    lines = [
        "# MMP9 DNA aptamer — synthesis shortlist", "",
        "## Stages completed", *[f"- {s}" for s in stage_note], "",
        "## Literature benchmarks (the bar)", bench, "",
        f"## Top {len(short)} candidates", "",
        "| # | sequence | len | tier | dG (kcal/mol) | consensus | jaccard | MFE | GC |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for i, r in enumerate(short, 1):
        lines.append(
            f"| {i} | `{r['sequence']}` | {r['length']} | {r['tier']} | "
            f"{r.get('dG_bind_kcal_mol','-')} | {r.get('consensus_score','-')} | "
            f"{r.get('contact_jaccard','-')} | {r['mfe']} | {r['gc_fraction']} |")
    lines += [
        "", "## How to read this",
        "- **tier 3** = survived all stages incl. MM/GBSA (most trustworthy);",
        "  tier 2 = consensus only; tier 1 = pilot oracle only.",
        "- **dG** is a RANKING signal (entropy neglected, implicit solvent) —",
        "  not a predicted Kd. Rank order is the usable output.",
        "- **jaccard** is contact-set agreement across independent predictions;",
        "  low values mean the pose is not reproducible even if scores look good.",
        "- **fnII_fraction** = share of contacts on the FnII exosite. For a pure",
        "  binder this is informational; for an inhibitor programme it is a",
        "  counter-selection criterion (nucleic acids there can ACTIVATE MMP9).",
        "", "## Next (wet lab)",
        "1. Synthesize the top candidates (5'-Cy5 label for MST).",
        "2. Measure Kd by MST or SPR; compare against F3B (~20 nM).",
        "3. Check selectivity against MMP2 and MMP7.",
        "4. Feed measured Kd back into the loop (active learning) for round 2.",
    ]
    (outdir / "REPORT.md").write_text("\n".join(lines) + "\n")

    print("\n".join(lines[:12]))
    print(f"\n-> {outdir}/shortlist.csv")
    print(f"-> {outdir}/REPORT.md")


if __name__ == "__main__":
    main()
