"""
MMP9 receptor-sequence provider for the Boltz-2 oracle (Stage 0 / handoff step 2).

The binding oracle co-folds the aptamer with the MMP9 protein, so it needs the
receptor amino-acid sequence — specifically the catalytic domain (UniProt P14780
residues 107-443), optionally with the FnII inserts removed to discourage the
activating FnII-exosite pose (see src/target/mmp9.py).

Integrity rule (same discipline as src/oracle/boltz2.py): this module NEVER
fabricates a protein sequence. It fetches the authoritative sequence from
UniProt when the network allows, or accepts one the user supplies; if neither is
available it raises with clear instructions. A wrong receptor sequence would
silently poison every downstream binding score, so guessing is not allowed.

Usage:
    # on an open-network workstation:
    seq = fetch_catalytic_domain()          # pulls P14780, slices 107-443
    # or, offline, from a FASTA the user pasted/downloaded:
    seq = catalytic_domain_from_fasta(open("P14780.fasta").read())
    # then:
    from oracle.boltz2 import Boltz2Oracle, Boltz2Config
    Boltz2Oracle(Boltz2Config(receptor_sequence=seq))
"""
from __future__ import annotations

import os
import ssl
import sys
import urllib.request
import urllib.error
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from target.mmp9 import UNIPROT, DOMAINS, LENGTH_AA  # noqa: E402

UNIPROT_FASTA_URL = f"https://rest.uniprot.org/uniprotkb/{UNIPROT}.fasta"
_CA_BUNDLE = os.environ.get("SSL_CERT_FILE") or "/root/.ccr/ca-bundle.crt"

_VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


def parse_fasta(text: str) -> str:
    """Return the sequence (concatenated, uppercase) from single-record FASTA."""
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    seq = "".join(ln for ln in lines if not ln.startswith(">")).upper()
    if not seq:
        raise ValueError("No sequence found in FASTA text.")
    bad = set(seq) - _VALID_AA
    if bad:
        raise ValueError(f"FASTA contains non-amino-acid characters: {sorted(bad)}")
    return seq


def _slice_domain(full_seq: str, span: tuple[int, int]) -> str:
    """Slice a 1-indexed inclusive UniProt residue range from the full sequence."""
    lo, hi = span
    if len(full_seq) < hi:
        raise ValueError(
            f"Sequence length {len(full_seq)} is shorter than the requested "
            f"range {span}; is this really the full P14780 preproprotein "
            f"(expected {LENGTH_AA} aa)?")
    return full_seq[lo - 1:hi]


def catalytic_domain_from_fasta(fasta_text: str, drop_fnii: bool = False) -> str:
    """Extract the catalytic domain (107-443) from a full-length P14780 FASTA.

    If drop_fnii=True, excise the FnII inserts (216-390) to make an FnII-free
    'mini' catalytic construct (mirrors the 2OVX experimental design), which
    discourages the activating FnII-exosite binding mode during co-folding.
    """
    full = parse_fasta(fasta_text)
    if len(full) != LENGTH_AA:
        # not fatal, but worth surfacing: isoforms / partial records differ
        print(f"[receptor] warning: fetched length {len(full)} != expected "
              f"{LENGTH_AA} for {UNIPROT}", file=sys.stderr)
    cat = _slice_domain(full, DOMAINS["catalytic_domain"])
    if not drop_fnii:
        return cat
    # remove the FnII span, re-based to the catalytic-domain start
    cat_start = DOMAINS["catalytic_domain"][0]
    f_lo, f_hi = DOMAINS["fnII_inserts"]
    i, j = f_lo - cat_start, f_hi - cat_start + 1
    return cat[:i] + cat[j:]


def fetch_full_sequence(timeout: int = 30) -> str:
    """Fetch the full P14780 sequence from UniProt. Raises (never fabricates) if
    the network/egress policy blocks it — the caller must then supply a FASTA."""
    ctx = ssl.create_default_context(
        cafile=_CA_BUNDLE if Path(_CA_BUNDLE).exists() else None)
    req = urllib.request.Request(UNIPROT_FASTA_URL,
                                 headers={"User-Agent": "findaptamer/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=timeout, context=ctx) as resp:
            return parse_fasta(resp.read().decode())
    except (urllib.error.URLError, urllib.error.HTTPError, ssl.SSLError, OSError) as e:
        raise RuntimeError(
            f"Could not fetch {UNIPROT} from UniProt ({e}). This session's egress "
            f"policy blocks {UNIPROT_FASTA_URL}. On an open-network workstation "
            f"this succeeds; otherwise download the FASTA and pass it to "
            f"catalytic_domain_from_fasta(). No sequence is fabricated.") from e


def fetch_catalytic_domain(drop_fnii: bool = False, timeout: int = 30) -> str:
    return catalytic_domain_from_fasta(_as_fasta(fetch_full_sequence(timeout)),
                                       drop_fnii=drop_fnii)


def _as_fasta(seq: str) -> str:
    return f">{UNIPROT}\n{seq}\n"


if __name__ == "__main__":
    print(f"Attempting to fetch MMP9 ({UNIPROT}) catalytic domain "
          f"({DOMAINS['catalytic_domain']})...")
    try:
        cat = fetch_catalytic_domain()
        print(f"catalytic domain ({len(cat)} aa):\n{cat}")
    except RuntimeError as e:
        print(f"\n[blocked] {e}\n")
        print("To proceed on the workstation:")
        print("  1. curl -o P14780.fasta "
              "'https://rest.uniprot.org/uniprotkb/P14780.fasta'")
        print("  2. python3 -c \"import sys; sys.path.insert(0,'src'); "
              "from target.receptor import catalytic_domain_from_fasta as f; "
              "print(f(open('P14780.fasta').read()))\"")
