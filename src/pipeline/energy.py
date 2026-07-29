"""
Stage 5 — physics energy tier: short MD + MM/GBSA binding-energy ranking.

This is the step AiDTA does not have at all: it ranks by docking/interface
heuristics and never computes an energy. MM/GBSA is the best-validated
quantitative ranker available for nucleic-acid complexes (moderate,
system-dependent correlation with experiment — see docs/02), so we use it to
convert "confident pose" into "ranked binding energy".

Method (single-trajectory MM/GBSA, the standard ranking protocol):
  1. Load the consensus complex structure, add hydrogens, parametrize with
     Amber14 (ff14SB protein + OL15 DNA) and GBn2 implicit solvent.
  2. Energy-minimize, then run a short MD at 300 K, collecting snapshots.
  3. For each snapshot compute
        dG_bind ~= E(complex) - E(protein) - E(DNA)
     with the SAME coordinates for all three (single-trajectory), which cancels
     most internal-energy noise. Solvation is included via the implicit model.
  4. Report the mean and standard error over snapshots.

Deliberate approximations (standard for ranking, stated openly):
  * configurational entropy (-TdS) is neglected — fine for RANKING congeneric
    binders, not for absolute affinities;
  * implicit solvent rather than explicit water + MMPBSA.py;
  * short sampling. Treat the output as a RANK, never as a predicted Kd.

Requires OpenMM on the GPU node:
    pip install openmm            (or: conda install -c conda-forge openmm)

Usage:
  python src/pipeline/energy.py --consensus results/consensus/consensus.json \
      --top 10 --ns 0.5 --out results/energy
"""
from __future__ import annotations
import argparse
import csv
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

NUCLEIC = {"DA", "DT", "DG", "DC", "A", "U", "G", "C", "T"}


def _require_openmm():
    try:
        import openmm  # noqa: F401
        from openmm import app  # noqa: F401
        import openmm.unit  # noqa: F401
    except ImportError as e:
        raise SystemExit(
            "OpenMM is required for the energy tier. On the GPU node run\n"
            "    pip install openmm\n"
            "or  conda install -c conda-forge openmm\n"
            f"({e})")


def split_complex(pdb_path: str, out_dir: Path):
    """Write protein-only and nucleic-only PDBs alongside the complex."""
    from openmm.app import PDBFile, Modeller

    pdb = PDBFile(str(pdb_path))
    top, pos = pdb.topology, pdb.positions

    def keep(which):
        m = Modeller(top, pos)
        drop = []
        for res in m.topology.residues():
            is_na = res.name.strip() in NUCLEIC or res.name.strip().rstrip("35") in NUCLEIC
            if (which == "protein" and is_na) or (which == "dna" and not is_na):
                drop.append(res)
        m.delete(drop)
        return m

    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {}
    for which in ("protein", "dna"):
        m = keep(which)
        p = out_dir / f"{which}.pdb"
        with open(p, "w") as f:
            PDBFile.writeFile(m.topology, m.positions, f)
        paths[which] = p
    return paths


def _system_from_pdb(pdb_path, forcefield, add_h=True):
    from openmm.app import PDBFile, Modeller
    pdb = PDBFile(str(pdb_path))
    modeller = Modeller(pdb.topology, pdb.positions)
    if add_h:
        modeller.addHydrogens(forcefield)
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=__import__("openmm.app", fromlist=["NoCutoff"]).NoCutoff,
        constraints=__import__("openmm.app", fromlist=["HBonds"]).HBonds,
    )
    return modeller, system


def potential_energy(modeller, system, positions=None):
    from openmm import LangevinMiddleIntegrator, Platform
    from openmm.app import Simulation
    import openmm.unit as unit
    integ = LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond,
                                     0.002 * unit.picoseconds)
    sim = Simulation(modeller.topology, system, integ)
    sim.context.setPositions(positions if positions is not None else modeller.positions)
    st = sim.context.getState(getEnergy=True)
    e = st.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
    del sim
    return e


def mmgbsa_single_structure(complex_pdb: str, work: Path, ns: float = 0.5,
                            snapshots: int = 20, minimize_steps: int = 2000):
    """Minimize + short implicit-solvent MD, then single-trajectory MM/GBSA."""
    from openmm import LangevinMiddleIntegrator
    from openmm.app import ForceField, Simulation, PDBFile, Modeller, NoCutoff, HBonds
    import openmm.unit as unit

    ff = ForceField("amber14-all.xml", "implicit/gbn2.xml")

    pdb = PDBFile(str(complex_pdb))
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(ff)
    system = ff.createSystem(modeller.topology, nonbondedMethod=NoCutoff,
                             constraints=HBonds)

    integ = LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond,
                                     0.002 * unit.picoseconds)
    sim = Simulation(modeller.topology, system, integ)
    sim.context.setPositions(modeller.positions)
    sim.minimizeEnergy(maxIterations=minimize_steps)

    total_steps = max(snapshots, int(ns * 1000 / 0.002))   # ns -> 2 fs steps
    stride = max(1, total_steps // snapshots)

    # identify protein vs nucleic atom indices for the decomposition
    na_idx, prot_idx = [], []
    for res in modeller.topology.residues():
        nm = res.name.strip().rstrip("35")
        tgt = na_idx if nm in NUCLEIC else prot_idx
        for a in res.atoms():
            tgt.append(a.index)
    if not na_idx or not prot_idx:
        raise RuntimeError(f"{complex_pdb}: could not find both protein and DNA")

    work.mkdir(parents=True, exist_ok=True)
    frames = []
    for i in range(snapshots):
        sim.step(stride)
        st = sim.context.getState(getPositions=True, getEnergy=True)
        frames.append(st.getPositions(asNumpy=True))

    # Build subsystems once, then evaluate each frame with the same coordinates.
    def subsystem(indices, name):
        m = Modeller(modeller.topology, modeller.positions)
        keep = set(indices)
        drop = [res for res in m.topology.residues()
                if not any(a.index in keep for a in res.atoms())]
        m.delete(drop)
        p = work / f"{name}.pdb"
        with open(p, "w") as f:
            PDBFile.writeFile(m.topology, m.positions, f)
        s = ff.createSystem(m.topology, nonbondedMethod=NoCutoff, constraints=HBonds)
        integ2 = LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond,
                                          0.002 * unit.picoseconds)
        return Simulation(m.topology, s, integ2), sorted(keep)

    sim_p, idx_p = subsystem(prot_idx, "protein")
    sim_d, idx_d = subsystem(na_idx, "dna")

    dGs = []
    for pos in frames:
        sim.context.setPositions(pos)
        e_c = sim.context.getState(getEnergy=True).getPotentialEnergy() \
            .value_in_unit(unit.kilocalorie_per_mole)
        sim_p.context.setPositions(pos[idx_p])
        e_p = sim_p.context.getState(getEnergy=True).getPotentialEnergy() \
            .value_in_unit(unit.kilocalorie_per_mole)
        sim_d.context.setPositions(pos[idx_d])
        e_d = sim_d.context.getState(getEnergy=True).getPotentialEnergy() \
            .value_in_unit(unit.kilocalorie_per_mole)
        dGs.append(e_c - e_p - e_d)

    n = len(dGs)
    mean = sum(dGs) / n
    var = sum((x - mean) ** 2 for x in dGs) / (n - 1) if n > 1 else 0.0
    return {"dG_bind_kcal_mol": round(mean, 2),
            "sem": round(math.sqrt(var / n), 2) if n > 1 else 0.0,
            "n_snapshots": n}


def pick_structure(rec: dict) -> str | None:
    """Best available predicted complex for a consensus record."""
    per = rec.get("per_predictor", {})
    best, best_score = None, -1.0
    for name, d in per.items():
        s = d.get("score", -1)
        st = d.get("structure")
        if st and Path(st).exists() and s > best_score:
            best, best_score = st, s
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--consensus", default="results/consensus/consensus.json")
    ap.add_argument("--top", type=int, default=10)
    ap.add_argument("--ns", type=float, default=0.5, help="MD length per candidate")
    ap.add_argument("--snapshots", type=int, default=20)
    ap.add_argument("--only-pass", action="store_true",
                    help="only candidates that passed consensus")
    ap.add_argument("--out", default="results/energy")
    args = ap.parse_args()

    _require_openmm()
    recs = json.loads(Path(args.consensus).read_text())
    if args.only_pass:
        recs = [r for r in recs if r.get("pass")]
    recs = recs[:args.top]
    if not recs:
        raise SystemExit("no candidates to process (try without --only-pass)")

    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    print(f"MM/GBSA on {len(recs)} candidates "
          f"({args.ns} ns, {args.snapshots} snapshots each)\n")

    rows = []
    for i, rec in enumerate(recs, 1):
        seq = rec["sequence"]
        struct = pick_structure(rec)
        if not struct:
            print(f"{i:3d}. {seq[:26]}..  no structure available — skipped")
            continue
        try:
            res = mmgbsa_single_structure(struct, outdir / f"cand{i:03d}",
                                          ns=args.ns, snapshots=args.snapshots)
        except Exception as e:
            print(f"{i:3d}. {seq[:26]}..  FAILED: {e}")
            continue
        row = {"sequence": seq, "length": len(seq),
               "consensus_score": rec.get("consensus_score"),
               "contact_jaccard": rec.get("contact_jaccard"), **res,
               "structure": struct}
        rows.append(row)
        print(f"{i:3d}. {seq[:26]}..  dG={res['dG_bind_kcal_mol']:>8.2f} "
              f"± {res['sem']:.2f} kcal/mol")

    rows.sort(key=lambda r: r["dG_bind_kcal_mol"])   # more negative = better
    (outdir / "energy.json").write_text(json.dumps(rows, indent=2, default=str))
    with open(outdir / "energy_ranked.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["rank", "sequence", "length", "dG_bind_kcal_mol", "sem",
                    "consensus_score", "contact_jaccard"])
        for k, r in enumerate(rows, 1):
            w.writerow([k, r["sequence"], r["length"], r["dG_bind_kcal_mol"],
                        r["sem"], r["consensus_score"], r["contact_jaccard"]])

    print(f"\nRanked by dG -> {outdir}/energy_ranked.csv")
    print("NOTE: MM/GBSA without entropy is a RANKING signal, not a predicted Kd.")


if __name__ == "__main__":
    main()
