"""Pseudogene-ORF specificity cut -- turn the non-novelty finding into a
quantified tumor-SPECIFICITY number, rigorously.

Context: 43/116 HCC pseudogene-ORF 'cryptic antigen' peptides are
exact substrings of the canonical proteome. Non-novelty is settled. The further
question: are they tumor-SPECIFIC? A peptide identical to a canonical-protein
substring is presented wherever that protein is -- so if the parent is expressed
in normal tissue, a T cell against it is on-target/off-tumor.

This script does it WITHOUT over-claiming:
  (1) identifies WHICH canonical protein(s) each peptide matches (parses GN= from
      the SwissProt header) -- so 'parent gene' is verified, not inferred from the
      locus name;
  (2) HARD evidence = direct observation in the HLA Ligand Atlas, an INDEPENDENT
      normal-tissue HLA-I ligandome in CANONICAL space (the correct reference for
      canonical-substring peptides -- IEAtlas-normal is cryptic-space and wrong
      here). Presence => presented on normal tissue => NOT tumor-specific. Absence
      is NOT proof of specificity (limited sampling) -> this is a conservative floor;
  (3) PRIOR (labeled, not measured) = parent is a ribosomal/housekeeping protein
      => ubiquitously expressed.

    python3 src/darkproteome/pseudogene_specificity.py

Public data only; stdlib + local SwissProt FASTA + HLA Ligand Atlas.
"""
import csv
import os
import re
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths  # noqa: E402  -- centralized data paths
import deepen_specificity as d  # SPROT, load_hla_ligand_atlas

csv.field_size_limit(10_000_000)

PRIMARY = os.path.join(paths.REPO, "data", "primary_tier1_nonnovelty.csv")
REAL = os.path.join(paths.REPO, "data", "claim_catalog_real.csv")
OUT = os.path.join(paths.REPO, "data", "pseudogene_specificity.csv")
HCC_DOI = "10.1126/sciadv.adn3628"
GN_RE = re.compile(r"\bGN=(\S+)")
PARENT_RE = re.compile(r"P\d+$")  # strip a trailing pseudogene index: RPS3AP12 -> RPS3A


def target_peptides():
    """the 43: HCC pseudogene-ORF peptides flagged canonical_self_exact=1."""
    pep = set()
    with open(PRIMARY, newline="") as fh:
        for r in csv.DictReader(fh):
            if (r["cohort"].startswith("HCC") and r["orf_class"] == "pseudogene-ORF"
                    and r["canonical_self_exact"] == "1"):
                pep.add(r["peptide"])
    return pep


def locus_map():
    """peptide -> pseudogene locus (orf_id_or_locus) from the HCC catalog rows."""
    m = {}
    with open(REAL, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r.get("citation_doi_pmid") == HCC_DOI and r.get("_canonical") == "no":
                s = (r.get("peptide_sequence") or "").strip().upper()
                if s and s.isalpha():
                    m.setdefault(s, r.get("orf_id_or_locus") or "?")
    return m


def matched_genes(peptides, sprot=d.SPROT):
    """peptide -> set of canonical gene symbols (GN=) of proteins that contain it.
    One streaming pass; for the ~43 short peptides a direct `in` test per protein
    is fast and exact."""
    hits = defaultdict(set)
    gene, seq = None, []

    def flush():
        if gene and seq:
            s = "".join(seq)
            for p in peptides:
                if p in s:
                    hits[p].add(gene)

    with open(sprot, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                mobj = GN_RE.search(line)
                gene = mobj.group(1) if mobj else line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())
        flush()
    return hits


def is_ribosomal(g):
    return bool(g) and (g.upper().startswith("RPS") or g.upper().startswith("RPL")
                        or g.upper().startswith("MRPL") or g.upper().startswith("MRPS"))


def main():
    paths.require(paths.SPROT, paths.HLALA)
    pep = target_peptides()
    loci = locus_map()
    print(f"target = HCC pseudogene-ORF canonical-self peptides: N={len(pep)}")

    print("identifying matched canonical gene(s) per peptide (GN= from SwissProt) ...")
    genes = matched_genes(pep)

    print("loading HLA Ligand Atlas (independent normal HLA-I ligandome, canonical space) ...")
    hla = set(d.load_hla_ligand_atlas())
    print(f"  HLA Ligand Atlas HLA-I peptides: {len(hla):,}")

    rows = []
    n_direct = n_ribo = n_parent_confirmed = n_any = 0
    for p in sorted(pep):
        gset = genes.get(p, set())
        locus = loci.get(p, "?")
        implied_parent = PARENT_RE.sub("", locus) if re.match(r"^[A-Z][A-Z0-9-]*P\d+$", locus or "") else ""
        parent_confirmed = bool(implied_parent) and implied_parent in gset
        ribo = any(is_ribosomal(g) for g in gset)
        in_hla = p in hla
        any_normal = in_hla or ribo  # hard OR strong-prior
        n_direct += in_hla
        n_ribo += ribo
        n_parent_confirmed += parent_confirmed
        n_any += any_normal
        rows.append((p, locus, ",".join(sorted(gset)) or "?", implied_parent or "-",
                     int(parent_confirmed), int(ribo), int(in_hla), int(any_normal)))

    n = len(pep)
    print(f"\n=== PSEUDOGENE SPECIFICITY (N={n}) ===")
    print(f"  matched canonical gene identified: {sum(1 for r in rows if r[2] != '?')}/{n}")
    print(f"  matched gene == named parent (verified, not inferred): {n_parent_confirmed}/{n}")
    print(f"  parent is RIBOSOMAL / housekeeping (strong normal-expression PRIOR): {n_ribo}/{n} "
          f"= {100*n_ribo/n:.0f}%")
    print(f"  DIRECTLY observed in HLA Ligand Atlas normal ligandome (HARD floor): {n_direct}/{n} "
          f"= {100*n_direct/n:.0f}%")
    print(f"  >>> ANY normal-presentation evidence (HLA-LA OR ribosomal parent): {n_any}/{n} "
          f"= {100*n_any/n:.0f}%")
    print(f"\n  SPECIFICITY NUMBER: of the {n} pseudogene 'cryptic antigen' peptides that are")
    print(f"  canonical substrings, {n_direct} are directly seen on normal tissue (hard floor) and")
    print(f"  {n_any} have normal-presentation evidence -> presumptively NOT tumor-specific.")
    print(f"  Of the original 116 HCC pseudogene claims that is {n_any}/116 = {100*n_any/116:.0f}%"
          f" flagged non-specific (lower bound; canonical-self was the gate).")

    # show the table
    print(f"\n  {'peptide':14s} {'locus':14s} {'matched_gene':16s} parent ribo hlaLA")
    for (p, locus, g, par, pc, ribo, hla_, anyn) in rows:
        print(f"  {p:14s} {locus:14s} {g:16s} {('Y' if pc else '-'):^6s} "
              f"{('Y' if ribo else '-'):^4s} {('Y' if hla_ else '-'):^5s}")

    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["peptide", "pseudogene_locus", "matched_canonical_genes", "implied_parent",
                    "parent_confirmed", "parent_ribosomal", "in_HLA_Ligand_Atlas_normal",
                    "any_normal_evidence"])
        w.writerows(rows)
    print(f"\nwrote -> {OUT}")
    print("\nCAVEATS: HLA-LA presence is a CONSERVATIVE floor (its sampling is finite; absence != "
          "specific). Ribosomal-parent is a labeled PRIOR, not measured expression -- see "
          "gtex_specificity.py for the measured GTEx parent-gene-expression extension "
          "covering all 43 parents.")


if __name__ == "__main__":
    main()
