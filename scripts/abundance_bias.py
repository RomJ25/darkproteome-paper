"""Test the abundance-bias hypothesis instead of asserting it -- and try hard to break it.

CORRECTION (ERROR #18). This docstring used to open by saying the catalogued
canonical-overlap rate (56.3%) "EXCEEDS" the library's latent ambiguity (nuORFdb, 34.1%), and treated
the difference as the thing to be explained. THAT COMPARISON IS INVALID AND IS WITHDRAWN. 56.3% is
over distinct catalogued PEPTIDES at native lengths, after search, FDR and dedup; 34.1% is over
distinct 9-MERS of an undetected candidate space. Different units, denominators and lengths -- the
same cross-unit error as the withdrawn "11-40x" fold-change. Do not quote a difference, ratio or
excess between those two numbers (nor against the 20.2% union or the 47.6% cap from library_union.py).

THE RESULTS BELOW ARE UNAFFECTED, because neither ever rested on that comparison. P1 is a
WITHIN-CATALOGUE contrast. P2 is a RATIO OF RATIOS -- an enrichment measured in the catalogue against
the same enrichment measured in the library -- which is dimensionless and so immune to the unit
mismatch. Only the framing sentence was wrong.

THE HYPOTHESIS, stated without the bad comparison: canonically-encoded peptides derive from abundant,
ubiquitously-expressed proteins, so they are over-detected in an immunopeptidome relative to their
share of the search space. The manuscript asserted this without testing it.

That is testable inside the resource. Two predictions, each with the control that could refute it.

P1  BREADTH OF DETECTION. A peptide of an abundant, ubiquitous protein should be detected across MORE
    of IEAtlas's 15 cancer types than a peptide with no canonical counterpart.
    CONFOUND: peptide length. Short peptides match the canonical proteome more readily by chance AND
    are the dominant, most-detectable HLA-I ligands -- length alone could produce the effect.
    CONTROL: stratify by length. The effect must hold WITHIN each length, or it is an artifact.

P2  ABUNDANT-CLASS ENRICHMENT. Ribosomal proteins (RPL*/RPS*) are the textbook abundant, ubiquitous
    housekeeping class. Detection bias predicts the canonical-overlapping epitopes are enriched for
    ORFs of ribosomal genes.
    CONFOUND -- AND THIS ONE IS FATAL IF UNCHECKED: the enrichment could be pure LIBRARY COMPOSITION.
    If nuORFdb's ribosomal-gene ORFs are simply more canonical-overlapping to begin with, we would
    see this enrichment with NO detection bias whatever, because it was in the search space already.
    CONTROL: measure the SAME enrichment in the library itself, over DISTINCT 9-mers -- the null world
    where nothing has been detected. Detection bias is the EXCESS of the catalogue's enrichment over
    the library's. If they agree, P2 measures composition, not detection, and we say so.

    CORRECTION (round-4 review). This used to measure the library side over raw (ORF,
    9-mer) OCCURRENCE PAIRS: every sliding window in every ORF counted separately, so a 9-mer repeated
    across N near-duplicate ORFs (ribosomal paralogs are exactly this) was counted N times. That is a
    DIFFERENT unit than the rest of the paper's library measurements (library_union.py's 34.1%/20.2%
    are over DISTINCT 9-mers, deduplicated globally) -- a reviewer pointed out that a ratio of ratios is
    dimensionless but not thereby unit-matched, and this was the concrete case of it. Recomputed over
    distinct 9-mers (a 9-mer counts once no matter how many ORFs contain it; "ribosomal" means it occurs
    in >=1 ribosomal-gene ORF): the library-side rate is 0.69x, not 0.91x -- STRONGER depletion, not
    weaker. The excess (2.51 / 0.69 = 3.65x) is therefore larger under the matched unit, not smaller.
    The qualitative conclusion survives the correction; the old 0.91x/2.77x numbers do not, and are
    retracted.

We can fail this. If the length-stratified effect vanishes, or the catalogue's ribosomal enrichment
merely reproduces the library's, the hypothesis does not survive and it comes OUT of the manuscript.

    python3 scripts/abundance_bias.py

NOTHING here bears on SOURCE. MS identifies the sequence, never the locus. Breadth of detection is
evidence about PRESENTATION, not about which locus produced the peptide.
"""
import csv
import json
import math
import os
import re
import statistics as st
import sys
from collections import defaultdict

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
EXT = os.path.join(REPO, "data", "external")
IE = os.path.join(EXT, "atlases", "IEAtlas_Epitopes_In_Cancer_Tissues.txt")
NUORF = os.path.join(EXT, "nuorfdb", "PA_nuORFdb_v1.2_protein.fasta")
SPROT = os.path.join(EXT, "swissprot_human.fasta")
TRANSLNC = os.path.join(EXT, "translnc", "lncRNA_peptide_AA_seq.fasta")
SCORED = os.path.join(REPO, "data", "claim_catalog_scored.csv")
csv.field_size_limit(10_000_000)
K = 9


def is_ribosomal(g):
    g = (g or "").upper()
    return g.startswith(("RPL", "RPS", "MRPL", "MRPS"))


_TRANSLNC_NONHUMAN_WGS_PREFIXES = ("AABR07", "CAAA01")
# Found by an independent, blind, from-scratch re-derivation of this whole computation (a
# round-5 check): the ALL-CAPS-means-human heuristic below also passes GenBank WGS-scaffold accessions,
# which are uppercase by GenBank convention regardless of species and are NOT HGNC gene symbols. Of the
# three distinct WGS-style prefixes present in this file, each was checked directly against its NCBI
# nuccore record: AABR07* (13,048 records) = Rattus norvegicus (rat), CAAA01* (3 records) = Mus musculus
# (mouse) -- both wrongly kept as "human" by the old heuristic -- but AUXG01* (53 records) is genuinely
# Homo sapiens (a human WGS scaffold, not a curated gene symbol, but still human) and must NOT be
# excluded. A blanket "exclude anything WGS-shaped" rule would itself be wrong here; only the two
# confirmed non-human prefixes are excluded.

def looks_human_translnc(header):
    """data/external/translnc/lncRNA_peptide_AA_seq.fasta is MULTI-SPECIES (found in a
    fresh-review pass) -- 435,173 human (ALL-CAPS gene-symbol convention, e.g. PAX8-AS1-...) headers
    followed by 148,667 mouse-convention headers (Tug1-..., Gm10619-..., RIKEN clone IDs). A header
    is 'human' iff the portion before its '-<digits>-<digits>aa' transcript/length marker contains no
    lowercase letters AND is not one of the confirmed-non-human WGS-scaffold prefixes (see
    _TRANSLNC_NONHUMAN_WGS_PREFIXES above). is_ribosomal(full_header) is provably equivalent to parsing
    the gene symbol first, since the gene symbol is always a strict prefix and is_ribosomal is a pure
    prefix test -- so no separate gene-token parser is needed for either function."""
    m = re.match(r"^(.*?)-\d+-\d+aa", header)
    prefix = m.group(1) if m else header
    if any(c.islower() for c in prefix):
        return False
    return not prefix.startswith(_TRANSLNC_NONHUMAN_WGS_PREFIXES)


def two_prop(k1, n1, k2, n2):
    p1, p2 = k1 / n1, k2 / n2
    pp = (k1 + k2) / (n1 + n2)
    se = math.sqrt(pp * (1 - pp) * (1 / n1 + 1 / n2)) or 1e-12
    return p1, p2, (p1 / p2 if p2 else float("inf")), (p1 - p2) / se


def ranksum_z(x, y):
    """Normal-approximation z for a Mann-Whitney rank-sum test; ties take midranks."""
    merged = sorted([(v, 0) for v in x] + [(v, 1) for v in y])
    ranks, i = {}, 0
    while i < len(merged):
        j = i
        while j + 1 < len(merged) and merged[j + 1][0] == merged[i][0]:
            j += 1
        r = (i + j) / 2 + 1
        for m in range(i, j + 1):
            ranks[m] = r
        i = j + 1
    n1, n2 = len(x), len(y)
    if not n1 or not n2:
        return 0.0
    u1 = sum(ranks[m] for m, (_, g) in enumerate(merged) if g == 0) - n1 * (n1 + 1) / 2
    sd = math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12) or 1e-12
    return (u1 - n1 * n2 / 2) / sd


def fasta(path):
    name, seq = None, []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].strip(), []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def main():
    for p in (IE, SCORED, NUORF, SPROT):
        if not os.path.exists(p):
            sys.exit(f"missing {p} (see data/SOURCES.md)")

    selfmap = {r["peptide"]: int(r["canonical_self"])
               for r in csv.DictReader(open(SCORED, newline="", encoding="utf-8"))}

    tissues, gene, length = defaultdict(set), {}, {}
    with open(IE, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.reader(fh, delimiter="\t")
        next(rd, None)
        for r in rd:
            if len(r) < 4 or not r[0]:
                continue
            p = r[0].strip().upper()
            if not p.isalpha():
                continue
            tissues[p].add((r[3] or "").strip())
            gene.setdefault(p, (r[2] or "").split("_")[0].upper())
            length.setdefault(p, len(p))

    scored = [p for p in tissues if p in selfmap]
    ov = [p for p in scored if selfmap[p]]
    nov = [p for p in scored if not selfmap[p]]
    n_types = len({t for s in tissues.values() for t in s})

    print(f"IEAtlas: {len(scored):,} scored cancer epitopes over {n_types} cancer types")
    print(f"  canonical-overlapping : {len(ov):,}")
    print(f"  NON-overlapping (ctrl): {len(nov):,}")

    # ---------------------------------------------------------------- P1
    print("\n" + "=" * 92)
    print("P1  breadth of detection -- are canonical-overlapping epitopes seen in MORE cancer types?")
    print("=" * 92)
    a = [len(tissues[p]) for p in ov]
    b = [len(tissues[p]) for p in nov]
    print(f"  {'':<27}{'n':>9}{'mean':>8}{'median':>8}{'>=2 types':>11}{'>=5 types':>11}")
    print("  " + "-" * 74)
    for lab, v in (("canonical-overlapping", a), ("NON-overlapping (control)", b)):
        g2 = 100 * sum(1 for x in v if x >= 2) / len(v)
        g5 = 100 * sum(1 for x in v if x >= 5) / len(v)
        print(f"  {lab:<27}{len(v):>9,}{st.mean(v):>8.2f}{st.median(v):>8.0f}{g2:>10.1f}%{g5:>10.1f}%")
    print(f"\n  pooled rank-sum z = {ranksum_z(a, b):,.0f}")

    print("\n  CONTROL -- stratified by peptide length (the confound: short peptides match canonical")
    print("  by chance more often, and are also the most detectable HLA-I ligands):\n")
    print(f"  {'len':>5}{'n(ovl)':>10}{'n(ctrl)':>10}{'mean(ovl)':>12}{'mean(ctrl)':>12}"
          f"{'z':>9}  holds?")
    print("  " + "-" * 74)
    holds, strata = 0, 0
    for L in sorted({length[p] for p in scored}):
        x = [len(tissues[p]) for p in ov if length[p] == L]
        y = [len(tissues[p]) for p in nov if length[p] == L]
        if len(x) < 30 or len(y) < 30:
            continue
        strata += 1
        z = ranksum_z(x, y)
        ok = st.mean(x) > st.mean(y)
        holds += ok
        print(f"  {L:>5}{len(x):>10,}{len(y):>10,}{st.mean(x):>12.2f}{st.mean(y):>12.2f}"
              f"{z:>9.1f}  {'yes' if ok else 'NO'}")
    p1_ok = holds == strata and strata >= 2
    print(f"\n  => effect holds in {holds}/{strata} length strata. "
          f"P1 {'SUPPORTED (not a length artifact)' if p1_ok else 'REFUTED / length-confounded'}.")

    # ---------------------------------------------------------------- P2
    print("\n" + "=" * 92)
    print("P2  enrichment for ribosomal-gene ORFs (the abundant housekeeping class)")
    print("=" * 92)
    ro = sum(1 for p in ov if is_ribosomal(gene.get(p)))
    rn = sum(1 for p in nov if is_ribosomal(gene.get(p)))
    c1, c2, c_rr, c_z = two_prop(ro, len(ov), rn, len(nov))
    print("  IN THE CATALOGUE (what was detected):")
    print(f"    canonical-overlapping      {ro:>8,}/{len(ov):<9,} = {100*c1:6.3f}%")
    print(f"    NON-overlapping (control)  {rn:>8,}/{len(nov):<9,} = {100*c2:6.3f}%")
    print(f"    risk ratio {c_rr:.2f}x   (z = {c_z:,.0f})")

    print("\n  THE CONTROL THAT COULD KILL IT -- the same enrichment in the LIBRARY, where nothing")
    print("  has been detected yet. Building canonical 9-mer set + scanning nuORFdb ...")
    canon = set()
    for _, s in fasta(SPROT):
        for i in range(len(s) - K + 1):
            canon.add(s[i:i + K])
    print(f"    canonical distinct {K}-mers: {len(canon):,}")

    # DISTINCT 9-mers, matching library_union.py's unit exactly -- NOT raw (ORF, 9-mer) occurrence
    # pairs. A 9-mer repeated across N near-duplicate ORFs (ribosomal paralogs are exactly this)
    # counts ONCE, same as a repeated catalogue peptide counts once. "Ribosomal" is set membership:
    # a 9-mer is ribosomal-associated iff it occurs in >=1 ribosomal-gene ORF (non-exclusive; a
    # kmer shared with a non-ribosomal ORF too is a small fraction -- see round-4 findings).
    ribo_kmers, nonribo_kmers = set(), set()
    n_orf = 0
    for h, s in fasta(NUORF):
        if len(s) < K:
            continue
        n_orf += 1
        m = re.search(r"GN=(\S+)", h)
        rib = is_ribosomal(m.group(1) if m else "")
        target = ribo_kmers if rib else nonribo_kmers
        for i in range(len(s) - K + 1):
            target.add(s[i:i + K])
    all_kmers = ribo_kmers | nonribo_kmers
    canon_kmers = all_kmers & canon
    noncanon_kmers = all_kmers - canon
    ov_ribo = sum(1 for k in canon_kmers if k in ribo_kmers)
    no_ribo = sum(1 for k in noncanon_kmers if k in ribo_kmers)
    l1, l2, l_rr, l_z = two_prop(ov_ribo, len(canon_kmers), no_ribo, len(noncanon_kmers))
    print(f"    scanned {n_orf:,} nuORFdb ORFs -> {len(all_kmers):,} distinct {K}-mers "
          f"({len(ribo_kmers):,} ribosomal-associated)")
    print(f"    canonical-overlapping      {ov_ribo:>8,}/{len(canon_kmers):<9,} = {100*l1:6.3f}%")
    print(f"    NON-overlapping            {no_ribo:>8,}/{len(noncanon_kmers):<9,} = {100*l2:6.3f}%")
    print(f"    risk ratio {l_rr:.2f}x   (z = {l_z:,.0f})   <-- LIBRARY-COMPOSITION baseline")

    excess = c_rr / l_rr if l_rr else float("inf")
    print(f"\n  catalogue {c_rr:.2f}x  vs  library {l_rr:.2f}x   ->  EXCESS = {excess:.2f}x")
    p2_ok = c_rr > l_rr * 1.2
    if p2_ok:
        print("  => P2 SUPPORTED: the catalogue is MORE ribosome-enriched than the search space it")
        print("     was drawn from. That excess cannot come from library composition; it arose during")
        print("     detection.")
    else:
        print("  => P2 REFUTED as evidence of DETECTION bias: the catalogue's enrichment is already")
        print("     present in the library. It measures composition, not detection. Report it as such.")

    # ---------------------------------------------------------------- Translnc-inclusive sensitivity
    # SCOPE, found in a fresh-review pass: the library baseline above is nuORFdb ALONE,
    # but IEAtlas's catalogue is drawn from all three of its integrated sources (nuORFdb + RPFdb +
    # Translnc). This was disclosed once in the Methods but never checked against the library this
    # paper already holds -- Translnc, whose union with nuORFdb was already computed for the headline
    # ambiguity rate (library_union.py: 34.1% -> 20.2%). Report both bounds here too, not just there.
    print("\n" + "=" * 84)
    print("TRANSLNC-INCLUSIVE SENSITIVITY (nuORFdb U Translnc, two variants)")
    print("=" * 84)
    translnc_variants = {}
    for key, human_only, label in (
        ("wholefile", False, "WHOLE Translnc file (matches library_union.py's own no-filter precedent)"),
        ("humanonly", True, "HUMAN-ONLY Translnc block (the more defensible choice for a human "
                             "immunopeptidome library)"),
    ):
        rk, nrk = set(ribo_kmers), set(nonribo_kmers)
        n_tl = 0
        for h, s in fasta(TRANSLNC):
            if len(s) < K:
                continue
            if human_only and not looks_human_translnc(h):
                continue
            n_tl += 1
            target = rk if is_ribosomal(h) else nrk
            for i in range(len(s) - K + 1):
                target.add(s[i:i + K])
        all_k = rk | nrk
        canon_k = all_k & canon
        noncanon_k = all_k - canon
        ov_r = sum(1 for k in canon_k if k in rk)
        no_r = sum(1 for k in noncanon_k if k in rk)
        tl1, tl2, tl_rr, tl_z = two_prop(ov_r, len(canon_k), no_r, len(noncanon_k))
        tl_excess = c_rr / tl_rr if tl_rr else float("inf")
        print(f"  {label} ({n_tl:,} peptides folded in):")
        print(f"    library RR {tl_rr:.2f}x   excess = {c_rr:.2f} / {tl_rr:.2f} = {tl_excess:.2f}x")
        translnc_variants[key] = {
            "n_translnc_peptides": n_tl, "library_rr": round(tl_rr, 2), "excess": round(tl_excess, 2),
        }
    print("\n  Both variants reported side by side, same as the 34.1%/20.2%/47.6% headline library")
    print("  numbers -- neither replaces the nuORFdb-only baseline above.")

    art = os.path.join(REPO, "data", "derived_detection_bias.json")
    json.dump({
        "n_cancer_types": n_types,
        "mean_types_overlapping": round(st.mean(a), 2),
        "mean_types_comparator": round(st.mean(b), 2),
        "pct_ge2_overlapping": round(100 * sum(1 for x in a if x >= 2) / len(a), 1),
        "pct_ge2_comparator": round(100 * sum(1 for x in b if x >= 2) / len(b), 1),
        "length_strata_tested": strata, "length_strata_holding": holds,
        "ribo_catalogue_rr": round(c_rr, 2),
        "ribo_library_rr": round(l_rr, 2),
        "ribo_excess": round(excess, 2),
        "ribo_pct_overlapping": round(100 * c1, 2),
        "ribo_pct_comparator": round(100 * c2, 2),
        "ribo_translnc_inclusive": translnc_variants,
    }, open(art, "w"), indent=2)
    print(f"\n  wrote {os.path.relpath(art, REPO)}")

    # ---------------------------------------------------------------- verdict
    print("\n" + "=" * 92)
    print("VERDICT")
    print("=" * 92)
    if p1_ok and p2_ok:
        print("""
  Both predictions survive their controls -- WITHOUT comparing the catalogue's 56.3% to the
  library's 34.1% as levels (ERROR #18; do not resurrect that comparison here). Canonically-encoded
  peptides are detected across more cancer types at every peptide length (P1, within-catalogue), and
  the catalogue's ribosomal enrichment (2.51x) is an excess OVER the library's own composition
  baseline (P2, a dimensionless ratio of ratios) -- BEYOND that class's share of the search space.

  THIS IS CORROBORATION, NOT PROOF, AND THE PAPER MUST SAY SO. Breadth of detection is a PROXY for
  abundance, not a measurement of it. A direct test needs protein-abundance data (PaxDb) and is not
  run here. The claim licensed is: "consistent with, and not explained by length or by library
  composition" -- nothing stronger.""")
    elif p1_ok:
        print("""
  P1 survives; P2 does not survive its control. Breadth of detection supports abundance bias, but the
  ribosomal enrichment is a property of the LIBRARY, not of detection -- it must be reported as
  composition, and must NOT be offered as evidence of detection bias.""")
    else:
        print("""
  P1 FAILS its control. The abundance-bias explanation is not supported by the resource's own data
  and must come OUT of the manuscript rather than stand as an untested hypothesis.""")
    print("""
  Scope, unchanged: MS identifies the SEQUENCE, never the LOCUS. This is evidence about what gets
  DETECTED, not about which locus produced any peptide.
""")
    return 0


if __name__ == "__main__":
    sys.exit(main())
