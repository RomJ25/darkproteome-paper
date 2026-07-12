"""What reference are we actually measuring against, and what would the ladder cost?

The 56.3% headline is a claim about SEQUENCE NOVELTY, and sequence novelty is REFERENCE-RELATIVE:
the estimand is N_i(R), never a bare N_i. The reviewer wants a nested sensitivity ladder --

    reviewed Swiss-Prot -> +reviewed isoforms -> +GENCODE coding translations
                        -> +common variants  -> +unreviewed UniProt

-- each layer with its exact release and hash, the exact and I/L-collapsed rates, and the
INCREMENTAL newly-matched count. The headline should sit on the NARROW, high-confidence layer and
be phrased "at least X%", because broader references can only ADD matches.

This script establishes exactly what is on disk, pins it, and states precisely what is missing.
It does not guess.

    python3 scripts/reference_provenance.py
"""
import hashlib
import os
import re
import sys

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
EXT = os.path.join(REPO, "data", "external")

# (layer, path-on-disk-or-None, what it is, why it matters)
LADDER = [
    ("1. reviewed Swiss-Prot (canonical)",
     os.path.join(EXT, "swissprot_human.fasta"),
     "reviewed human canonical proteins, one sequence per gene",
     "HIGH-CONFIDENCE. The conservative floor: a match here is unambiguous."),
    ("2. + reviewed isoforms",
     None,
     "UniProt human reviewed + isoform sequences (uniprot_sprot_varsplic)",
     "An ncORF peptide matching an ISOFORM of a canonical gene is not novel either."),
    ("3. + GENCODE protein-coding translations",
     None,
     "GENCODE pc_translations.fa (the annotated coding space)",
     "Broader annotated coding space; catches genes Swiss-Prot has not reviewed."),
    ("4. + common variants",
     None,
     "e.g. gnomAD-derived variant peptides",
     "A peptide differing by one common SNP is not a novel ORF product."),
    ("5. + unreviewed UniProt (TrEMBL)",
     None,
     "the full predicted/unreviewed protein universe",
     "MOST INCLUSIVE. Detects more competing sources -- but is no longer cleanly a "
     "'canonical protein' rate, so it CANNOT carry the headline."),
]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fasta_stats(path):
    n, res = 0, 0
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
            else:
                res += len(line.strip())
    return n, res


def main():
    print("=" * 78)
    print("REFERENCE PROVENANCE — the estimand is N_i(R), never a bare N_i")
    print("=" * 78)

    have, missing = [], []
    for name, path, what, why in LADDER:
        if path and os.path.exists(path):
            n, res = fasta_stats(path)
            print(f"\n  [ON DISK] {name}")
            print(f"    file    : {os.path.relpath(path, REPO)}")
            print(f"    sha256  : {sha256(path)}")
            print(f"    bytes   : {os.path.getsize(path):,}")
            print(f"    entries : {n:,} sequences, {res:,} residues")
            print(f"    is      : {what}")
            print(f"    role    : {why}")
            have.append(name)
        else:
            print(f"\n  [ MISSING ] {name}")
            print(f"    would be: {what}")
            print(f"    role    : {why}")
            missing.append(name)

    print("\n" + "=" * 78)
    print("WHAT WE CAN AND CANNOT SAY TODAY")
    print("=" * 78)
    print(f"  layers on disk : {len(have)} of {len(LADDER)}")
    print(f"  layers missing : {len(missing)} of {len(LADDER)}")
    print()
    print("  The published 56.3% is measured against LAYER 1 ONLY. That is the conservative,")
    print("  high-confidence layer, so the honest phrasing is already available:")
    print()
    print('      "at least 56.3% of catalogued cryptic cancer epitopes in IEAtlas are exact')
    print('       substrings of reviewed canonical human proteins (Swiss-Prot, release pinned)."')
    print()
    print("  'AT LEAST' is load-bearing: broader references can only ADD matches, never remove")
    print("  them, so every missing layer can push the rate UP and none can push it down. The")
    print("  headline is therefore SAFE as a lower bound even with the ladder unbuilt --")
    print("  but the SENSITIVITY CLAIM ('and it does not depend on the reference') is NOT yet")
    print("  supported, and must not be made until layers 2-5 are measured.")
    print()
    print("  What the missing layers change:")
    print("    - they can only INCREASE the overlap rate;")
    print("    - so they STRENGTHEN 'the record cannot distinguish the ncORF from a canonical")
    print("      source by sequence' (the paper's actual claim);")
    print("    - and they WEAKEN any statement about the size of a 'genuinely novel remainder'")
    print("      (which is why that phrasing is retired).")
    print()
    print("  => The direction of the unmeasured risk is KNOWN and runs in the paper's favour.")
    print("     That is worth stating explicitly rather than leaving a referee to wonder.")

    if missing:
        print("\n  TO BUILD THE LADDER, these must be downloaded (all public, all documented in")
        print("  data/SOURCES.md conventions):")
        for m in missing:
            print(f"    - {m}")


if __name__ == "__main__":
    main()
