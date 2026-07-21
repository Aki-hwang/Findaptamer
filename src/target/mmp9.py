"""
MMP9 target definition for structure-guided aptamer design.

Human matrix metalloproteinase-9 (gelatinase B), UniProt P14780, 707-aa
preproenzyme. This module encodes the domain map, the druggable/epitope
definitions, the mandatory FnII counter-selection zone, the recommended
docking receptors, and the known-aptamer benchmarks — everything the
generator/oracle stages need to stay on-target.

Sources (compiled 2026-07): UniProt P14780; Rowsell 2002 J Mol Biol 319:173
(1GKC); Elkins 2002 (1L6J); Tochowicz 2007 J Mol Biol 371:989 (2OVX series);
Appleby 2017 J Biol Chem PMID 28235803 (5TH6 apo, 5TH9 andecaliximab Fab);
Shimada 2018 Biochem J 475:1597 (nucleic-acid potentiation via FnII);
De Franciscis/Cerchia WO2013153138A1 + PLOS ONE 2016 (F3B RNA aptamer);
LVMH/INSERM US9902961B2 (G-quadruplex inhibitory DNA aptamers).

CRITICAL BIOLOGY: single-stranded DNA/RNA binds the MMP9 fibronectin type-II
(FnII) exosite (~res 216-390) and *potentiates* (activates) catalysis. An
inhibitory aptamer MUST be steered away from FnII and toward the catalytic
cleft or the validated allosteric (andecaliximab) patch. FNII_COUNTERSELECT
below is used to penalize/counter-select FnII-docking candidates.
"""
from __future__ import annotations

UNIPROT = "P14780"
LENGTH_AA = 707

# Domain map (UniProt numbering, full preproprotein). Ranges are approximate
# where the literature disagrees on exact splits (flagged).
DOMAINS = {
    "signal_peptide":      (1, 19),      # cleaved
    "prodomain":           (20, 106),    # cysteine-switch; Cys99 ligates Zn in zymogen
    "catalytic_domain":    (107, 443),   # interrupted by FnII inserts
    "fnII_inserts":        (216, 390),   # gelatin-binding exosite (approx module splits)
    "og_linker":           (444, 511),   # O-glycosylated flexible hinge
    "hemopexin_PEX":       (512, 707),   # four-blade beta-propeller
}

MATURE_N_TERMINUS = 107  # Phe107 after activation cleavage Arg106-Phe107

# Catalytic machinery
CATALYTIC = {
    "zinc_ligands": [401, 405, 411],   # His401, His405, His411 (HEXXHXXGXXH)
    "catalytic_base": 402,             # Glu402
    "cysteine_switch": 99,             # Cys99 (zymogen latency)
    "s1_prime_pocket": [421, 423],     # Pro421-X-Tyr423 (selectivity determinant)
    "specificity_loop": (424, 430),    # Arg424-Pro430
}

# ---- Candidate epitopes an aptamer can be designed against --------------------
# Each epitope lists key surface residues (UniProt numbering) that define the
# "desired binding site", mirroring how AiDTA defines an epitope from an
# antibody-antigen interface.
EPITOPES = {
    # Clinically validated, selective, NON-COMPETITIVE inhibitory site.
    # Defined from the andecaliximab (GS-5745) Fab-MMP9 complex (PDB 5TH9).
    # Best size/mechanism fit for a ~30-40 nt aptamer; away from FnII.
    "allosteric_andecaliximab": {
        "residues": [162, 111, 113, 198],   # R162 (most critical), E111, D113, I198
        "mechanism": "non-competitive allosteric inhibition + blocks pro-MMP9 activation",
        "defined_from_pdb": "5TH9",
        "recommended_receptor": "5TH6",      # apo MMP9 (no ligand bias)
        "goal": "inhibitor",
        "notes": "Near a Ca2+ pocket at the prodomain-catalytic junction.",
    },
    # Competitive active-site / S1' pocket. Mechanistically direct but the cleft
    # is deep/narrow -> harder for a bulky oligonucleotide to occlude.
    "active_site_s1prime": {
        "residues": [401, 405, 411, 402, 421, 423, 424, 425, 426, 427, 428, 429, 430],
        "mechanism": "competitive active-site / S1' occlusion",
        "defined_from_pdb": "1GKC",
        "recommended_receptor": "1GKC",      # 2.3 A, wild-type, S1' defined by inhibitor
        "goal": "inhibitor",
        "notes": "Use 2OVX-series (mini domain, no FnII) for cleanest S1' SAR.",
    },
}

# HARD CONSTRAINT: counter-select against binding here (activation risk).
FNII_COUNTERSELECT = {
    "residues_range": DOMAINS["fnII_inserts"],
    "reason": "Nucleic acids binding the FnII exosite POTENTIATE MMP9 (Shimada 2018).",
    "action": "penalize/counter-select candidates whose predicted pose contacts FnII",
}

# Recommended docking receptors (fetch coordinates from RCSB when proxy allows).
RECEPTORS = {
    "5TH6": "apo human MMP9 (no ligand bias) — allosteric/exosite docking",
    "1L6J": "proMMP9 N-terminal (prodomain+catalytic+FnII) — prodomain-junction targeting",
    "1GKC": "MMP9 catalytic domain + FnII + reverse-hydroxamate inhibitor, 2.3 A, WT — active site",
    "2OVX": "mini catalytic domain (FnII removed), E402Q, high-res — clean S1' pocket",
}

# ---- Known aptamer benchmarks (the bar to beat) ------------------------------
# affinity_nM: reported Kd; function: what it actually does.
BENCHMARK_APTAMERS = [
    {
        "name": "F3B",
        "type": "RNA (2'-F pyrimidines; F3Bomf variant 2'-OMe purines)",
        "sequence": "UGCCCUGCCCUCACCCGUUAGCCUGAGCGCCCCGCA",
        "length": 36,
        "affinity_nM": 20.0,
        "selective_vs_MMP2_7": True,
        "function": "binder/imaging (inhibition NOT established)",
        "source": "WO2013153138A1; PLOS ONE 2016 10.1371/journal.pone.0149387",
    },
    {
        "name": "A3",
        "type": "ssDNA",
        "sequence": None,   # 40-mer; exact sequence not in accessible sources
        "length": 40,
        "affinity_nM": 13.4,
        "selective_vs_MMP2_7": None,
        "function": "ACTIVATOR (potentiates catalysis) — a cautionary benchmark, not a goal",
        "source": "Shimada 2018 Biochem J 475:1597 (PMC5941315)",
    },
    {
        "name": "LVMH_G4_series",
        "type": "ssDNA G-quadruplex",
        "sequence": None,   # multiple SEQ IDs in patent; exact winner needs PDF
        "length": None,
        "affinity_nM": None,
        "selective_vs_MMP2_7": None,
        "function": "INHIBITOR (gelatinase inhibition claimed) — closest prior art to our goal",
        "source": "US9902961B2 / JP6599874B2 (LVMH Recherche & INSERM)",
    },
]


def summary() -> dict:
    return {
        "uniprot": UNIPROT,
        "length_aa": LENGTH_AA,
        "domains": DOMAINS,
        "catalytic": CATALYTIC,
        "epitopes": list(EPITOPES.keys()),
        "fnII_counterselect": FNII_COUNTERSELECT["residues_range"],
        "receptors": list(RECEPTORS.keys()),
        "benchmark_best_nM": min(a["affinity_nM"] for a in BENCHMARK_APTAMERS
                                 if a["affinity_nM"]),
    }


if __name__ == "__main__":
    import json
    print(json.dumps(summary(), indent=2, default=str))
