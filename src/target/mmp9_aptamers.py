"""
Curated, provenance-tracked set of published MMP-9 aptamers.

Why this file exists: an in-silico SELEX needs a fitness signal, and every
fitness signal we can build has to be *calibrated against real binders* before
it is allowed to rank anything. `src/target/mmp9.py` carried three benchmark
entries but only one had a usable sequence (F3B, a 2'-F/2'-OMe RNA), which made
every calibration run a test of modified-RNA modelling as much as of binding.

Compiling this set fixed that and turned up a G-quadruplex prior: the DNA-side
selections against MMP-9 give G4 folds. That is the strongest target-specific
prior we have, and unlike a predicted binding energy it is an observation
rather than a model output.

PROVENANCE CAVEAT, and it is load-bearing. F3B, 8F14A and the LVMH patent all
come from ONE laboratory -- Toulme and Dausse at INSERM Bordeaux appear as
authors or inventors on all three. So the G4 observations are NOT independent
selections; they share a lab, and plausibly a library and a protocol. A shared
protocol bias is a live alternative explanation for the convergence. Weigh the
prior accordingly, and see the `independent_of` field on each entry.

CONFIDENCE field — how much weight a calibration may put on an entry:
  "high"    sequence verified in >=2 independent primary sources
  "medium"  sequence from one primary source
  "claimed" property asserted in the literature but the sequence itself is not
            in any source reachable from here (do NOT use as a positive control)

Compiled 2026-08 by web search; every URL below was retrieved, not inferred.
Where a value could not be retrieved it is None rather than a guess.
"""
from __future__ import annotations

DNA = "DNA"
RNA_2F = "RNA (2'-F pyrimidine)"

KNOWN_APTAMERS = [
    {
        "name": "F3B",
        "chemistry": RNA_2F,
        "sequence": "UGCCCUGCCCUCACCCGUUAGCCUGAGCGCCCCGCA",
        "length": 36,
        "affinity_nM": 20.0,
        "fold": "stem-loop (no G4 reported)",
        "forms_g4": False,
        "function": "binder / tumour imaging (inhibition NOT established)",
        "confidence": "medium",
        "lab": "Toulme / Dausse (INSERM Bordeaux)",
        "lab_group": "toulme_dausse",
        "independent_of": [],   # same lab as 8F14A and LVMH_G4_series
        "notes": (
            "Truncated to 36 nt; purines substituted with 2'-O-methyl in the "
            "nuclease-resistant F3Bomf variant. Selected with a 2'-F-pyrimidine "
            "RNA library against recombinant hMMP-9. Kd by filter retention. "
            "BINDING DOMAIN ON MMP-9 IS NOT MAPPED in any source we could reach "
            "-- this is the confounder that makes it a weak positive control "
            "for a catalytic-domain receptor (see docs/06)."
        ),
        "sources": [
            "WO2013153138A1 (De Franciscis / Cerchia / Toulme)",
            "PLOS ONE 2016 10.1371/journal.pone.0149387",
        ],
    },
    {
        "name": "8F14A",
        "chemistry": RNA_2F,
        "sequence": None,          # not in any source reachable from here
        "length": None,
        "affinity_nM": None,
        "fold": "G-quartet (confirmed by CD and Tm)",
        "forms_g4": True,           # experimental (CD + Tm)
        "function": "binder at a SECOND epitope, non-competitive with F3B",
        "confidence": "claimed",
        "lab": "Toulme / Dausse (INSERM Bordeaux)",
        "lab_group": "toulme_dausse",
        "independent_of": [],   # same lab as F3B and LVMH_G4_series
        "notes": (
            "Used as the sandwich partner to F3B in a piezoelectric biosensor: "
            "the two aptamers bind MMP-9 simultaneously without competition, so "
            "MMP-9 presents at least two independent aptamer epitopes. "
            "G-quartet formation confirmed experimentally -> evidence point 2 "
            "for the G4 convergence."
        ),
        "sources": [
            "Scarano, Dausse, Crispo, Toulme, Minunni, "
            "Anal. Chim. Acta 2015, 897, 1-9",
        ],
    },
    {
        "name": "MMP9-DNA-30",
        "chemistry": DNA,
        "sequence": "TCGTATGGCACGGGGTTGGTGTTGGGTTGG",
        "length": 30,
        "affinity_nM": None,       # Kd not stated in the sources we reached
        "fold": "G-rich; G4 predicted (see src/analysis/g4.py)",
        "forms_g4": True,           # COMPUTATIONAL prediction only
        "function": "binder / biosensor capture element",
        "confidence": "high",
        "lab": None,            # original selection not traced to a source we could reach
        "lab_group": None,      # unknown -- counted as its own group, optimistically
        "independent_of": ["F3B", "8F14A", "LVMH_G4_series"],
        "notes": (
            "THE MOST USEFUL ENTRY: unmodified DNA, so it can be modelled and "
            "synthesised exactly as published -- no 2'-F / 2'-OMe confounder. "
            "Re-used unchanged by at least three independent groups across "
            "photoacoustic, transistor and SERS platforms over 2021-2026, "
            "which is strong empirical evidence that it binds. "
            "Deployed with a 5'-phosphate-T10 spacer in the IGZO work; the "
            "spacer is immobilisation chemistry, not part of the aptamer. "
            "CAUTION: its 3' half is close to the thrombin-binding aptamer "
            "(see g4.py) -- quantify that before treating it as MMP-9-specific."
        ),
        "sources": [
            "Baik et al., Photoacoustics 2021, PMC8521288 "
            "(doi:10.1016/j.pacs.2021.100294)",
            "ACS Appl. Mater. Interfaces 2025, 10.1021/acsami.5c03926 "
            "(tear MMP-9, IGZO TFT)",
            "Anal. Bioanal. Chem. 2026, 10.1007/s00216-026-06573-4 "
            "(SERS plasmonic ruler)",
        ],
    },
    {
        "name": "A3",
        "chemistry": "ssDNA",
        "sequence": None,
        "length": 40,
        "affinity_nM": 13.4,
        "fold": None,
        "forms_g4": None,           # not reported
        "function": "ACTIVATOR -- potentiates catalysis via the FnII exosite",
        "lab": "Shimada et al.",
        "lab_group": "shimada",
        "independent_of": ["F3B", "8F14A", "LVMH_G4_series"],
        "confidence": "claimed",
        "notes": (
            "A cautionary benchmark, not a goal: tight binding to MMP-9 by "
            "ssDNA can ACTIVATE the enzyme. Any designed sequence must be "
            "counter-selected against the FnII exosite "
            "(mmp9.FNII_COUNTERSELECT)."
        ),
        "sources": ["Shimada 2018 Biochem J 475:1597 (PMC5941315)"],
    },
    {
        "name": "LVMH_G4_series",
        "chemistry": DNA,
        "sequence": None,          # patent SEQ IDs not reachable from here
        "length": None,
        "affinity_nM": None,
        "fold": "G-quadruplex (explicitly claimed in the patent)",
        "forms_g4": True,           # asserted in the patent
        "function": "INHIBITOR of MMP-9 gelatinase activity; cell-penetrant",
        "lab": "Toulme / Dausse (INSERM Bordeaux) + LVMH Recherche",
        "lab_group": "toulme_dausse",
        "independent_of": [],   # same lab as F3B and 8F14A
        "inventors": ["Eric Dausse", "Jean-Jacques Toulme",
                      "Jean Hubert Cauchard", "Robin Kurfurst",
                      "Sylvianne Schnebert"],
        "confidence": "claimed",
        "notes": (
            "Selected from a phosphotriester oligonucleotide library for "
            "sequences that both bind MMP-9 and abolish its enzymatic "
            "activity. The patent states the winning aptamer has a "
            "G-quadruplex structure -> evidence point 3 for G4 convergence. "
            "Sequences are behind patent-office and Google-Patents endpoints "
            "that refuse requests from this sandbox; retrievable from a normal "
            "browser."
        ),
        "sources": [
            "US9902961B2 / US20160326530A1 / FR3015986A1 "
            "(LVMH Recherche & INSERM)",
        ],
    },
]


# ---- reference G-quadruplex aptamers, for the specificity sanity check -------
# These are NOT MMP-9 binders. They are the two most famous G4 aptamers, and
# they are here because G4 aptamers are known to recur across unrelated SELEX
# campaigns. If a "novel MMP-9" design is really just one of these, we need to
# know before anyone synthesises it.
REFERENCE_G4_APTAMERS = [
    {
        "name": "TBA / HD1 / ARC-183",
        "target": "thrombin exosite I",
        "sequence": "GGTTGGTGTGGTTGG",
        "fold": "antiparallel chair, two G-tetrads, TGT + 2x TT loops",
        "sources": ["Bock et al. Nature 1992; PDB 148D/4DII"],
    },
    {
        "name": "AS1411 / AGRO100",
        "target": "nucleolin",
        "sequence": "GGTGGTGGTGGTTGTGGTGGTGGTGG",
        "fold": "parallel/mixed G4, polymorphic",
        "sources": ["Bates et al. Exp Mol Pathol 2009"],
    },
]


def with_sequence(chemistry: str | None = None, min_confidence: str = "medium"):
    """Entries usable as positive controls: sequence known, confidence high enough.

    `min_confidence` "medium" admits medium+high; "high" admits high only.
    "claimed" entries are never returned -- they have no sequence to score.
    """
    order = {"claimed": 0, "medium": 1, "high": 2}
    floor = order[min_confidence]
    return [a for a in KNOWN_APTAMERS
            if a["sequence"]
            and order[a["confidence"]] >= floor
            and (chemistry is None or a["chemistry"] == chemistry)]


def g4_evidence():
    """Observations behind the G-quadruplex prior, WITH their independence.

    Reported as a dict rather than a list because the count alone overstates
    the case: two of the three G4 observations come from the same laboratory,
    so "three observations" is not "three independent selections".
    """
    # Match on an explicit flag, never on the prose of `fold`: F3B's fold reads
    # "stem-loop (no G4 reported)", so a substring test for "G4" counted a
    # NEGATIVE observation as supporting evidence.
    hits = [a for a in KNOWN_APTAMERS if a.get("forms_g4") is True]
    # Count canonical group ids, not the free-text lab names: those differ by a
    # trailing "+ LVMH Recherche" and would be tallied as two groups.
    labs = {a["lab_group"] for a in hits if a["lab_group"]}
    return {
        "observations": [
            {"aptamer": a["name"], "chemistry": a["chemistry"],
             "evidence": a["fold"], "confidence": a["confidence"],
             "lab": a["lab"], "source": a["sources"][0]}
            for a in hits
        ],
        "n_observations": len(hits),
        "n_distinct_labs": len(labs) + sum(1 for a in hits if not a["lab"]),
        "caveat": (
            "Toulme/Dausse (INSERM Bordeaux) authored or invented both 8F14A "
            "and the LVMH series, so those two share a lab and plausibly a "
            "library and protocol. The only G4 observation outside that group "
            "is MMP9-DNA-30, and that one is a computational prediction "
            "(G4Hunter), not an experimental structure. Shared protocol bias "
            "remains a live alternative explanation for the convergence."
        ),
    }


if __name__ == "__main__":
    import json
    print(json.dumps({
        "n_entries": len(KNOWN_APTAMERS),
        "usable_as_positive_control": [a["name"] for a in with_sequence()],
        "dna_positive_controls": [a["name"] for a in with_sequence(DNA)],
        "g4_evidence": g4_evidence(),
    }, indent=2))
