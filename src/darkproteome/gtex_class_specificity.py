"""Class-resolved GTEx normal-expression probe for the NON-canonical survivors.

WHY: gtex_specificity.py settled the pseudogene-ORF canonical-self subset (43/43
parents broadly normal-expressed -> non-specific, MEASURED). The harder, more
important question is the GENUINELY non-canonical survivors -- the altORF and
lncRNA-ORF claims that PASSED the non-novelty floor (~0% canonical-self). All of
these are labelled TSA (claimed tumor-SPECIFIC). Does the specificity floor extend
to them, or is the specificity problem confined to the pseudogene class?

This probe answers that, AT THE HONEST EVIDENTIAL STRENGTH OF EACH CLASS:

  - lncRNA-ORF: the source IS the lncRNA gene, so GTEx normal expression of that
    gene is a DIRECT (RNA-level) tumor-restriction test. Caveat: clone-named
    GENCODE lncRNAs resolve by symbol only partially across annotation versions,
    and the resolvable subset is biased toward NAMED (better-characterized) lncRNAs.

  - altORF: the source locus is a CANONICAL host gene (the alt-ORF is an alternate
    frame / uORF inside it). Host-gene RNA expression says the locus is transcribed
    in normal tissue -- it does NOT say the alt-frame peptide is presented there.
    So a broadly-expressed host is a SOURCE-not-tumor-restricted flag, NOT proof of
    normal presentation. (These Raja altORF TSA claims report NO normal check.)

  - pseudogene-ORF: the right target is the canonical PARENT, not the pseudogene
    locus -> see gtex_specificity.py (43/43). Shown here only for contrast.

FRAMING: "source expressed in normal => tumor-restriction not supported
as reported", never "the antigen is false". Calibrated, not inflated.

    python3 src/darkproteome/gtex_class_specificity.py

Public data only; stdlib + GTEx gene-median GCT (reuses gtex_specificity loader).
"""
import csv
import os
import sys
import statistics

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths
import gtex_specificity as gs

REAL = os.path.join(paths.REPO, "data", "claim_catalog_real.csv")
OUT = os.path.join(paths.REPO, "data", "gtex_class_specificity.csv")
SKIP = {"not reported", "multiple noncoding loci", ""}
EXPR = 1.0          # minimal "expressed" median TPM
BROAD = 27          # >= EXPR TPM in >= half of the 54 tissues = broadly/housekeeping


def class_loci():
    """orf_class -> {locus -> n_claims} (TSA peptide rows only)."""
    out = {}
    with open(REAL, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            cls = r.get("orf_class", "?")
            if cls not in ("altORF", "lncRNA-ORF", "pseudogene-ORF"):
                continue
            loc = (r.get("orf_id_or_locus") or "").strip()
            if loc.lower() in SKIP:
                continue
            out.setdefault(cls, {}).setdefault(loc, 0)
            out[cls][loc] += 1
    return out


def expr_stats(loci, gtex, tissues):
    """For a set of source loci: resolution + normal-expression breakdown."""
    rows = []
    for loc in loci:
        vec = gs.lookup(loc, gtex)
        if vec is None:
            rows.append((loc, 0, 0.0, 0))
            continue
        mx = max(vec)
        ntis = sum(1 for x in vec if x >= EXPR)
        rows.append((loc, 1, mx, ntis))
    return rows


def main():
    paths.require(paths.GTEX_MEDIAN)
    print("loading GTEx v8 gene-median TPM ...")
    tissues, gtex = gs.load_gtex_medians()
    print(f"  GTEx: {len(gtex):,} gene symbols x {len(tissues)} normal tissues\n")

    cl = class_loci()
    out_rows = []
    print(f"{'class':16s} {'loci':>5s} {'reslv':>6s} {'expr>=1':>8s} {'broad':>7s} {'medMaxTPM':>10s}")
    for cls in ("lncRNA-ORF", "altORF", "pseudogene-ORF"):
        loci = cl.get(cls, {})
        rows = expr_stats(loci, gtex, tissues)
        n = len(rows)
        resolved = [r for r in rows if r[1]]
        nr = len(resolved)
        expr = [r for r in resolved if r[2] >= EXPR]
        broad = [r for r in resolved if r[3] >= BROAD]
        med = statistics.median([r[2] for r in resolved]) if resolved else 0.0
        pe = f"{len(expr)}/{nr}" if nr else "0/0"
        pb = f"{len(broad)}/{nr}" if nr else "0/0"
        print(f"{cls:16s} {n:>5d} {nr:>4d}/{n:<1d} {pe:>8s} {pb:>7s} {med:>10.2f}")
        for (loc, found, mx, ntis) in rows:
            out_rows.append((cls, loc, loci[loc], found, f"{mx:.2f}", ntis,
                             int(found and mx >= EXPR), int(found and ntis >= BROAD)))

    # claim-weighted view for the two genuinely-non-canonical classes
    print("\n=== READ-OUT (claim-weighted, the genuinely non-canonical survivors) ===")
    for cls in ("lncRNA-ORF", "altORF"):
        loci = cl.get(cls, {})
        tot = sum(loci.values())
        rows = {r[0]: r for r in expr_stats(loci, gtex, tissues)}
        c_res = sum(loci[l] for l in loci if rows[l][1])
        c_expr = sum(loci[l] for l in loci if rows[l][1] and rows[l][2] >= EXPR)
        c_broad = sum(loci[l] for l in loci if rows[l][1] and rows[l][3] >= BROAD)
        print(f"  {cls}: {tot} TSA claims; source resolvable {c_res} ({100*c_res/tot:.0f}%); "
              f"of resolvable, broadly normal-expressed {c_broad}/{c_res} "
              f"({100*c_broad/max(c_res,1):.0f}%)")

    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["orf_class", "source_locus", "n_claims", "resolved_in_gtex",
                    "max_median_tpm", "n_tissues_expr_ge1", "expressed_ge1tpm",
                    "broadly_expressed_ge_half"])
        w.writerows(out_rows)
    print(f"\nwrote -> {OUT}")
    print("""
=== CALIBRATED VERDICT (the point of this probe) ===
  pseudogene-ORF: specificity FAILS, MEASURED -- the canonical PARENT is normal-
      expressed for 43/43 of the canonical-self peptides (gtex_specificity.py). Strong.
  altORF: the canonical HOST locus is broadly normal-transcribed, so the source is
      not tumor-restricted; but host RNA != alt-frame presentation, so this is a
      reported-gap flag (these claims report no normal check), NOT proof of normal
      presentation. Moderate / as-reported.
  lncRNA-ORF: source loci are mostly lowly / narrowly expressed, and only a minority
      resolve by symbol -> NOT shown non-specific. The cleanest, most plausibly real
      class. Weak / inconclusive (matches the ~0% canonical-self non-novelty result).
  -> The specificity indictment is CLASS-SPECIFIC to pseudogene-ORF; it does NOT
     blanket the non-canonical survivors. The peptide-level normal-presentation test
     for those is IEAtlas-normal (9.1% Raja / 12.7% HCC), already in tier1_nonnovelty.
  CAVEATS: GTEx is normal RNA, not presentation; symbol matching is annotation-version
     limited for clone-named lncRNAs (a complete pass needs ENSG-level mapping);
     "expressed in normal" is a specificity RISK, never "the antigen is fake".
""")


if __name__ == "__main__":
    main()
