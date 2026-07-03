"""GTEx parent-gene expression extension to the pseudogene-ORF specificity floor.

Turns the LABELED ribosomal-family PRIOR (pseudogene_specificity.py) into a MEASURED
normal-tissue expression number, using GTEx v8 gene-median TPM (56,200 genes x 54
normal tissues).

Each of the 43 HCC
pseudogene-ORF 'cryptic antigen' peptides is an EXACT substring of one or more
canonical proteins (matched_canonical_genes, verified via SwissProt GN=). A peptide
identical to a canonical-protein substring is presented wherever that protein is
made. So if a matched canonical parent gene is EXPRESSED IN NORMAL TISSUE, a T cell
against the peptide is on-target / off-tumor -> the peptide is presumptively NOT
tumor-specific. GTEx gives the normal-tissue expression directly.

This is the MEASURED complement to the two evidence tiers already in
pseudogene_specificity.py:
  - HLA Ligand Atlas direct observation  = HARD floor (normal PRESENTATION), and
  - ribosomal-parent family               = a labeled PRIOR (assumed expression).
GTEx replaces that prior with measured RNA expression across 54 normal tissues, and
covers ALL 43 parents, not just the ribosomal ones.

CONSERVATIVE BY CONSTRUCTION (no over-claiming):
  - normal-tissue medians (not tumor-adjacent); we report several TPM thresholds
    (>=1 standard "expressed"; >=10 robustly expressed) plus tissue BREADTH;
  - "expressed in normal" => specificity RISK, never "the peptide is fake";
  - absence from GTEx is NOT proof of specificity (symbol aliasing, ncORF parents);
  - when a peptide matches several canonical genes we take the MAX expression (any
    expressed parent already makes the peptide non-specific).

    python3 src/darkproteome/gtex_specificity.py

DATA SOURCE (re-download if data/external/gtex/ is ever lost; ~7 MB, open access):
  GTEx v8 gene-median TPM (genes x 54 tissues), GTEx Analysis 2017-06-05:
  https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz
  -> data/external/gtex/GTEx_v8_gene_median_tpm.gct.gz   (resolved by paths.GTEX_MEDIAN)

Public data only; stdlib + local GTEx gene-median GCT (data/external/gtex/),
reusing the verified peptide->parent-gene mapping from pseudogene_specificity.py.
"""
import csv
import gzip
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths
import pseudogene_specificity as ps
import deepen_specificity as d

OUT = os.path.join(paths.REPO, "data", "gtex_pseudogene_specificity.csv")
EXPR = 1.0    # standard minimal "expressed" median TPM in a tissue
ROBUST = 10.0  # robustly expressed

# A few SwissProt GN= symbols differ from the GENCODE v26 symbols GTEx ships.
# Map the handful that occur among the 43 parents to their GTEx-v26 alias so the
# measurement is not lost to a naming change (confirmed 1:1 via HGNC).
ALIAS = {
    "STMP1": "C7orf73", "ARK2N": "C18orf25", "LLPH": "C12orf32",
    "ATP5PB": "ATP5F1",  # HGNC ATP5PB == GENCODE-v26 ATP5F1 (ENSG00000116459)
    "H2BC1": "HIST1H2BA", "H2BC3": "HIST1H2BB", "H2BC4": "HIST1H2BC",
    "H2BC5": "HIST1H2BD", "H2BC6": "HIST1H2BE", "H2BC7": "HIST1H2BF",
    "H2BC8": "HIST1H2BG", "H2BC9": "HIST1H2BH", "H2BC10": "HIST1H2BI",
    "H2BC11": "HIST1H2BJ", "H2BC12": "HIST1H2BK", "H2BC13": "HIST1H2BL",
    "H2BC14": "HIST1H2BM", "H2BC15": "HIST1H2BN", "H2BC17": "HIST1H2BO",
    "H2BC18": "HIST1H2BPS1", "H2BC21": "HIST2H2BE", "H2BC12L": "HIST3H2BB",
}


def load_gtex_medians(path=paths.GTEX_MEDIAN):
    """gene symbol (UPPER) -> [median TPM per tissue]; elementwise max if a symbol
    appears on more than one ENSG row. Returns (tissue_names, sym2vec)."""
    sym2vec = {}
    with gzip.open(path, "rt") as fh:
        fh.readline()                                # '#1.2'
        fh.readline()                                # 'nrows<TAB>ncols'
        tissues = fh.readline().rstrip("\n").split("\t")[2:]
        for line in fh:
            f = line.rstrip("\n").split("\t")
            sym = f[1].strip().upper()
            vec = [float(x) for x in f[2:]]
            if sym in sym2vec:
                sym2vec[sym] = [max(a, b) for a, b in zip(sym2vec[sym], vec)]
            else:
                sym2vec[sym] = vec
    return tissues, sym2vec


def lookup(gene, gtex):
    """GTEx vector for a gene symbol, trying the HGNC<->GENCODE-v26 alias."""
    g = gene.upper()
    if g in gtex:
        return gtex[g]
    a = ALIAS.get(gene, "").upper()
    return gtex.get(a)


def main():
    paths.require(paths.GTEX_MEDIAN, paths.SPROT, paths.HLALA)

    pep = ps.target_peptides()
    loci = ps.locus_map()
    print(f"target = HCC pseudogene-ORF canonical-self peptides: N={len(pep)}")
    print("identifying matched canonical gene(s) per peptide (SwissProt GN=) ...")
    genes = ps.matched_genes(pep)
    print("loading HLA Ligand Atlas (normal HLA-I ligandome, the hard floor) ...")
    hla = set(d.load_hla_ligand_atlas())
    print("loading GTEx v8 gene-median TPM ...")
    tissues, gtex = load_gtex_medians()
    print(f"  GTEx: {len(gtex):,} gene symbols x {len(tissues)} normal tissues")

    rows = []
    n = len(pep)
    n_found = n_expr = n_robust = 0
    n_old = n_new = 0           # HLA-LA OR prior  vs  HLA-LA OR GTEx-expressed
    n_hla = 0
    breadth = []                # # tissues expressed, per peptide (for the median)
    for p in sorted(pep):
        gset = genes.get(p, set())
        vecs = [v for g in gset if (v := lookup(g, gtex)) is not None]
        found = bool(vecs)
        if found:
            per_tissue = [max(v[i] for v in vecs) for i in range(len(tissues))]
            max_tpm = max(per_tissue)
            n_tis = sum(1 for x in per_tissue if x >= EXPR)
        else:
            max_tpm, n_tis = 0.0, 0
        expr = found and max_tpm >= EXPR
        robust = found and max_tpm >= ROBUST
        in_hla = p in hla
        ribo = any(ps.is_ribosomal(g) for g in gset)
        old_any = in_hla or ribo
        new_any = in_hla or expr
        n_found += found
        n_expr += expr
        n_robust += robust
        n_hla += in_hla
        n_old += old_any
        n_new += new_any
        if expr:
            breadth.append(n_tis)
        rows.append((p, loci.get(p, "?"), ",".join(sorted(gset)) or "?",
                     int(found), f"{max_tpm:.1f}", n_tis, int(in_hla),
                     int(ribo), int(expr), int(robust), int(new_any)))

    med_breadth = sorted(breadth)[len(breadth) // 2] if breadth else 0
    print(f"\n=== GTEx PARENT-GENE EXPRESSION (N={n}) ===")
    print(f"  matched parent gene resolvable in GTEx:                 {n_found}/{n}")
    print(f"  parent expressed in normal (>= {EXPR:.0f} TPM, >=1 tissue):     {n_expr}/{n} "
          f"= {100*n_expr/n:.0f}%")
    print(f"  parent ROBUSTLY expressed (>= {ROBUST:.0f} TPM):                {n_robust}/{n} "
          f"= {100*n_robust/n:.0f}%")
    print(f"  median # of 54 normal tissues the parent is expressed in: {med_breadth} "
          f"(of the expressed ones) -> housekeeping breadth")
    print(f"\n  --- effect on the specificity floor ---")
    print(f"  any-normal, HLA-LA hard OR ribosomal PRIOR (labeled):    {n_old}/{n} "
          f"= {100*n_old/n:.0f}%   (= {n_old}/116 = {100*n_old/116:.0f}% of HCC pseudogene claims)")
    print(f"  any-normal, HLA-LA hard OR GTEx MEASURED >= {EXPR:.0f} TPM:        {n_new}/{n} "
          f"= {100*n_new/n:.0f}%   (= {n_new}/116 = {100*n_new/116:.0f}% of HCC pseudogene claims)")
    print(f"  (HLA-LA hard floor alone: {n_hla}/{n} = {100*n_hla/n:.0f}%)")

    print(f"\n  {'peptide':18s} {'matched_gene':20s} {'maxTPM':>7s} {'nTis':>4s} hla ribo expr")
    for (p, locus, g, found, mtpm, ntis, in_hla, ribo, expr, robust, newany) in rows:
        print(f"  {p:18s} {g[:20]:20s} {mtpm:>7s} {ntis:>4d}  "
              f"{('Y' if in_hla else '-'):^3s} {('Y' if ribo else '-'):^4s} "
              f"{('Y' if expr else '-'):^4s}")

    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["peptide", "pseudogene_locus", "matched_canonical_genes",
                    "gene_in_gtex", "max_median_tpm", "n_tissues_expr_ge1",
                    "in_HLA_Ligand_Atlas_normal", "parent_ribosomal",
                    "gtex_expressed_ge1tpm", "gtex_robust_ge10tpm",
                    "any_normal_evidence_measured"])
        w.writerows(rows)
    print(f"\nwrote -> {OUT}")
    print("\nCAVEATS: GTEx medians are normal-tissue RNA, not protein presentation; >=1 TPM is the "
          "standard minimal-expression call (we also report >=10 and tissue breadth). 'Expressed in "
          "normal' => specificity RISK (on-target/off-tumor), NOT 'the peptide is false'. Absence "
          "from GTEx is not proof of specificity. The HLA Ligand Atlas tier remains the HARD floor; "
          "GTEx supersedes the ribosomal PRIOR with a measured number over all 43 parents.")


if __name__ == "__main__":
    main()
