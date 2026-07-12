"""Test the abundance-bias hypothesis instead of asserting it -- and try hard to break it.

The catalogued canonical-overlap rate (56.3%) EXCEEDS the latent ambiguity of the library IEAtlas
searched (nuORFdb, 34.1%). The manuscript states, but does NOT test, an explanation: detection bias.
Canonically-encoded peptides derive from abundant, ubiquitously-expressed proteins, so they are
over-detected in an immunopeptidome relative to their share of the search space.

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
    CONTROL: measure the SAME enrichment in the library itself, over (ORF, 9-mer) candidate-epitope
    pairs -- the null world where nothing has been detected. Detection bias is the EXCESS of the
    catalogue's enrichment over the library's. If they agree, P2 measures composition, not detection,
    and we say so.

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
SCORED = os.path.join(REPO, "data", "claim_catalog_scored.csv")
csv.field_size_limit(10_000_000)
K = 9


def is_ribosomal(g):
    g = (g or "").upper()
    return g.startswith(("RPL", "RPS", "MRPL", "MRPS"))


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

    lib = {(True, True): 0, (True, False): 0, (False, True): 0, (False, False): 0}
    n_orf = 0
    for h, s in fasta(NUORF):
        if len(s) < K:
            continue
        n_orf += 1
        m = re.search(r"GN=(\S+)", h)
        rib = is_ribosomal(m.group(1) if m else "")
        for i in range(len(s) - K + 1):
            lib[(s[i:i + K] in canon, rib)] += 1
    n_ov = lib[(True, True)] + lib[(True, False)]
    n_no = lib[(False, True)] + lib[(False, False)]
    l1, l2, l_rr, l_z = two_prop(lib[(True, True)], n_ov, lib[(False, True)], n_no)
    print(f"    scanned {n_orf:,} nuORFdb ORFs -> {n_ov + n_no:,} (ORF, {K}-mer) candidate epitopes")
    print(f"    canonical-overlapping      {lib[(True, True)]:>8,}/{n_ov:<9,} = {100*l1:6.3f}%")
    print(f"    NON-overlapping            {lib[(False, True)]:>8,}/{n_no:<9,} = {100*l2:6.3f}%")
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
    }, open(art, "w"), indent=2)
    print(f"\n  wrote {os.path.relpath(art, REPO)}")

    # ---------------------------------------------------------------- verdict
    print("\n" + "=" * 92)
    print("VERDICT")
    print("=" * 92)
    if p1_ok and p2_ok:
        print("""
  Both predictions survive their controls. The catalogued overlap rate (56.3%) exceeding the
  library's latent ambiguity (34.1%) is CONSISTENT WITH detection bias: canonically-encoded peptides
  are detected across more cancer types at every peptide length, and the catalogue is enriched for
  the abundant housekeeping class BEYOND that class's share of the search space.

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
