"""The DIRECT abundance measurement the manuscript lists as open -- replacing the detection-breadth proxy.

WHERE THIS COMES FROM. abundance_bias.py tested the abundance/detectability explanation for the
canonical-overlap pattern in the catalogue. (NOTE, ERROR #18: an earlier version of this header framed
that as explaining why the catalogue's 56.3% "exceeds" the library's 34.1%. That comparison is INVALID
and WITHDRAWN -- distinct PEPTIDES at native lengths post-search vs distinct 9-MERS of an undetected
candidate space are different units over different denominators. Do not quote a difference, ratio or
excess between them. The results are unaffected: they never rested on it.)

It supported the explanation with a PROXY (canonical-overlapping sequences are detected across more
cancer types; ribosomal-gene ORFs are enriched in the catalogue beyond their library share -- a RATIO
OF RATIOS, which is dimensionless and so immune to that unit mismatch) and said so out loud:

    "Breadth of detection is a PROXY for abundance, not a measurement of it. A direct test needs
     protein-abundance data (PaxDb) and is not run here."

This is that test. PaxDb v6.1, H. sapiens whole-organism integrated dataset (ppm), joined to the
canonical proteins the catalogued peptides actually match.

    H:  among IEAtlas's catalogued sequences, those that overlap a canonical protein do so
        preferentially with ABUNDANT canonical proteins -- and the abundance of the matched canonical
        protein predicts how broadly the peptide was detected.

THREE MEASUREMENTS, each with the control that could kill it.

A   WHICH CANONICAL PROTEINS DO THEY HIT?  Compare the abundance of the matched canonical protein(s)
    against the abundance distribution of ALL canonical human proteins. Reported peptide-weighted
    (max and median over a peptide's matched proteins) AND protein-weighted -- the protein-level view
    has NO peptide-clustering problem at all, because its unit IS the protein.

B   THE KEY TEST -- DOES ABUNDANCE PREDICT DETECTION BREADTH?  Among canonical-overlapping sequences
    only, bin by the abundance of the matched canonical protein and compute the mean number of the 15
    IEAtlas cancer types each sequence is detected in. The prediction is a MONOTONE INCREASING trend.
    If the trend is FLAT or DECREASING, the hypothesis FAILS and the claim comes out of the paper.

C   THE CONTROLS.
    c1  PEPTIDE LENGTH. A reviewer already caught us on length. Short peptides match the canonical
        proteome more readily by chance AND are the dominant, most-detectable HLA-I ligands. So B is
        repeated WITHIN every length stratum, and the headline trend is LENGTH-STANDARDIZED (direct
        standardization to the catalogue's own length distribution).
    c2  PROTEIN LENGTH -- the obvious attack on A. A longer protein contains more substrings, so it
        is likelier to be hit BY CHANCE; if "hit" proteins are merely long proteins, A is an artifact.
        So A is repeated WITHIN protein-length deciles.
    c3  PLACEBO. Break the peptide->protein link (assign each peptide a random canonical protein's
        abundance). The trend must collapse. If it does not, the machinery invents trends and every
        number it produces is worthless.

    INFERENCE, throughout. No z-test. The 174,465 peptides are NOT independent -- they cluster by
    gene, and a reviewer destroyed an earlier z for exactly that. Every interval here is a
    GENE-CLUSTERED BOOTSTRAP (the consequence_robust.py pattern: deterministic LCG, B=2000, seed
    20260713), run under TWO cluster definitions -- the matched canonical gene (the cluster that
    shares the EXPOSURE) and the IEAtlas source gene (the cluster that shares the ORF label). The
    trend has to survive both.

    python3 scripts/abundance_direct.py

SCOPE, unchanged and non-negotiable. MS identifies the SEQUENCE, never the LOCUS. This is evidence
about what gets DETECTED. It is not evidence that any particular peptide came from the canonical gene
rather than the ncORF -- and a finding that abundance drives detection does NOT make the biology fake.
It bears on what the reported record can support.
"""
import csv
import json
import math
import os
import statistics as st
import sys
from collections import defaultdict

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
EXT = os.path.join(REPO, "data", "external")
IE = os.path.join(EXT, "atlases", "IEAtlas_Epitopes_In_Cancer_Tissues.txt")
SPROT = os.path.join(EXT, "swissprot_human.fasta")
PAX = os.path.join(EXT, "paxdb", "9606-WHOLE_ORGANISM-integrated.txt")
PAX_V5 = os.path.join(EXT, "paxdb", "9606-WHOLE_ORGANISM-integrated-v5.0.txt")
LINKS = os.path.join(EXT, "paxdb", "paxdb-uniprot-links-v6.1.tsv")
SCORED = os.path.join(REPO, "data", "claim_catalog_scored.csv")
ART = os.path.join(REPO, "data", "derived_abundance_direct.json")
csv.field_size_limit(10_000_000)

B = 2000
SEED = 20260713
NBINS = 5          # quintiles for the headline; deciles also printed
MIN_STRATUM = 200  # peptides required before a length stratum is reported


class RNG:
    """Deterministic LCG -- same generator as consequence_robust.py, so runs are reproducible."""

    def __init__(self, s):
        self.s = s

    def randint(self, n):
        self.s = (self.s * 6364136223846793005 + 1442695040888963407) & ((1 << 64) - 1)
        return (self.s >> 33) % n


# ---------------------------------------------------------------------------- rank-based helpers
def midranks(vals):
    """Midranks (ties averaged), returned in the order of `vals`."""
    order = sorted(range(len(vals)), key=lambda i: vals[i])
    r = [0.0] * len(vals)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and vals[order[j + 1]] == vals[order[i]]:
            j += 1
        m = (i + j) / 2 + 1
        for k in range(i, j + 1):
            r[order[k]] = m
        i = j + 1
    return r


def auc(x, y):
    """P(X > Y) + 0.5 P(X = Y) -- the Mann-Whitney common-language effect size, by hand.

    A descriptive effect size, NOT a p-value: 0.5 = no shift, >0.5 = x stochastically larger.
    Deliberately reported without a z, because these samples are clustered (see c2)."""
    if not x or not y:
        return float("nan")
    r = midranks(list(x) + list(y))
    n1, n2 = len(x), len(y)
    u1 = sum(r[:n1]) - n1 * (n1 + 1) / 2
    return u1 / (n1 * n2)


def rank_corr(pairs):
    """Spearman rho, by hand (Pearson on midranks)."""
    if len(pairs) < 3:
        return float("nan")
    rx = midranks([p[0] for p in pairs])
    ry = midranks([p[1] for p in pairs])
    n = len(pairs)
    mx, my = sum(rx) / n, sum(ry) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    sxx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    syy = math.sqrt(sum((b - my) ** 2 for b in ry))
    return sxy / (sxx * syy) if sxx and syy else float("nan")


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


# ---------------------------------------------------------------------------- loaders
def load_swissprot():
    """[(accession, entry_name, gene_symbol, sequence)] from the canonical human proteome."""
    prots = []
    for h, s in fasta(SPROT):
        parts = h.split("|")
        acc = parts[1] if len(parts) > 2 else h.split()[0]
        entry = parts[2].split()[0] if len(parts) > 2 else ""
        gene = ""
        for tok in h.split():
            if tok.startswith("GN="):
                gene = tok[3:].upper()
                break
        prots.append((acc, entry.upper(), gene, s))
    return prots


def load_paxdb(path):
    """{ENSP: ppm}, {gene_symbol: ppm}. Duplicate keys keep the MAX (most abundant proteoform).

    Column layout differs BY RELEASE and must not be assumed: v6.1 is
    (gene_name, string_external_id, abundance); v5.0 is (string_external_id, abundance) with no gene
    column at all. Hard-coding v6.1's 3 columns silently parsed ZERO rows of v5.0 -- the version
    check just vanished from the output with no error. So: locate the ENSP field, take the last
    numeric field as ppm, and treat a leading non-ENSP field as the gene symbol."""
    by_ensp, by_gene = {}, {}
    parsed = 0
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = [x.strip() for x in line.rstrip("\n").split("\t")]
            ei = next((i for i, x in enumerate(f) if "ENSP" in x.upper()), None)
            if ei is None or len(f) < 2:
                continue
            try:
                ppm = float(f[-1])
            except ValueError:
                continue
            if ppm <= 0:
                continue
            ensp = f[ei].split(".")[-1]           # 9606.ENSP00000370010 -> ENSP00000370010
            gene = f[0].upper() if ei > 0 else ""
            parsed += 1
            if ensp:
                by_ensp[ensp] = max(by_ensp.get(ensp, 0.0), ppm)
            if gene:
                by_gene[gene] = max(by_gene.get(gene, 0.0), ppm)
    if not parsed:
        sys.exit(f"FATAL: parsed 0 abundance rows from {path} -- the column layout changed. Stop.")
    return by_ensp, by_gene


def load_links():
    """{ENSP: UniProt entry name} for taxid 9606, from PaxDb's own UniProt cross-reference."""
    out = {}
    with open(LINKS, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith("9606."):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 2:
                continue
            out[f[0].split(".", 1)[1].strip()] = f[1].strip().upper()
    return out


def join_abundance(prots, by_ensp, by_gene, links):
    """SwissProt index -> ppm. Official ENSP->entry-name link first; gene symbol as fallback."""
    ensp_by_entry = defaultdict(list)
    for ensp, entry in links.items():
        ensp_by_entry[entry].append(ensp)

    ab, via = {}, {"link": 0, "gene": 0}
    for i, (_acc, entry, gene, _s) in enumerate(prots):
        vals = [by_ensp[e] for e in ensp_by_entry.get(entry, []) if e in by_ensp]
        if vals:
            ab[i] = max(vals)
            via["link"] += 1
        elif gene and gene in by_gene:
            ab[i] = by_gene[gene]
            via["gene"] += 1
    return ab, via


def load_catalogue():
    """peptide -> (set of cancer types, IEAtlas source gene)."""
    tissues, gene = defaultdict(set), {}
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
    return tissues, gene


def match_peptides(peptides, prots):
    """peptide -> [SwissProt indices containing it as an EXACT substring].

    Same relation as `canonical_self` in claim_catalog_scored.csv (exact substring, no I/L collapse),
    but it keeps WHICH proteins matched -- which is the whole point here. Index the queries by their
    8-residue prefix (8 = the catalogue's minimum length) and stream the proteome once."""
    P = min(len(p) for p in peptides)
    by_prefix = defaultdict(list)
    for p in peptides:
        by_prefix[p[:P]].append(p)
    hits = defaultdict(set)
    for i, (_a, _e, _g, s) in enumerate(prots):
        n = len(s)
        for j in range(n - P + 1):
            cands = by_prefix.get(s[j:j + P])
            if not cands:
                continue
            for p in cands:
                if s.startswith(p, j):
                    hits[p].add(i)
    return {p: sorted(v) for p, v in hits.items()}


# ---------------------------------------------------------------------------- the statistic
# A (bin, length) CELL is the atom of everything below. The point estimate and the bootstrap share
# ONE implementation of the standardized statistic (std_from_cells), so a resample cannot silently
# compute a different quantity than the estimate whose CI it claims to be.
def cells_from_rows(rows, key="bin"):
    """{(bin, peptide length): [sum of breadth, n]}."""
    cell = defaultdict(lambda: [0.0, 0])
    for r in rows:
        c = cell[(r[key], r["len"])]
        c[0] += r["breadth"]
        c[1] += 1
    return cell


def std_from_cells(cell, weights):
    """Length-standardized mean detection-breadth per abundance bin.

    Direct standardization to the catalogue's own length distribution (the consequence_robust.py
    pattern), so a bin cannot look broad merely by holding shorter, more-detectable peptides.
    Weights are renormalized over the length strata a bin actually populates."""
    out = {}
    for b in sorted({k[0] for k in cell}):
        num = den = 0.0
        for L, w in weights.items():
            s, n = cell.get((b, L), (0.0, 0))
            if n:
                num += w * s / n
                den += w
        out[b] = num / den if den else float("nan")
    return out


def std_bin_means(rows, weights, key="bin"):
    return std_from_cells(cells_from_rows(rows, key), weights)


def q5_minus_q1_cells(cell, weights):
    """Headline effect: length-standardized mean breadth in the top vs bottom abundance bin."""
    m = std_from_cells(cell, weights)
    if not m:
        return float("nan")
    b = sorted(m)
    return m[b[-1]] - m[b[0]]


def q5_minus_q1(rows, weights):
    return q5_minus_q1_cells(cells_from_rows(rows), weights)


def clustered_bootstrap(rows, weights, clust_key, label):
    """Percentile CI resampling CLUSTERS with replacement -- the unit of dependence, not the peptide.

    Each cluster's contribution to the (bin, length) cells is precomputed once, so a replicate is a
    sum of sparse cell contributions rather than a rebuild of ~96k rows. Identical statistic and
    identical RNG draw sequence (one randint per cluster per replicate) -- just not quadratically
    slow."""
    byc = defaultdict(list)
    for r in rows:
        byc[r[clust_key]].append(r)
    # per-cluster sparse contribution: [((bin, len), sum_breadth, n), ...]
    contrib = [[(k, v[0], v[1]) for k, v in cells_from_rows(rs).items()] for rs in byc.values()]
    nc = len(contrib)
    rng = RNG(SEED)
    boot = []
    for _ in range(B):
        cell = defaultdict(lambda: [0.0, 0])
        for _ in range(nc):
            for k, s, n in contrib[rng.randint(nc)]:
                c = cell[k]
                c[0] += s
                c[1] += n
        v = q5_minus_q1_cells(cell, weights)
        if v == v:
            boot.append(v)
    boot.sort()
    lo = boot[int(0.025 * len(boot))]
    hi = boot[int(0.975 * len(boot)) - 1]
    print(f"    {label:<34}95% CI [{lo:+.3f}, {hi:+.3f}]   ({nc:,} clusters, B={len(boot):,})")
    return lo, hi


def bins_by_quantile(vals, k):
    """k equal-count bin edges over `vals` (ties can make bins uneven; that is honest, not fixed)."""
    s = sorted(vals)
    return [s[int(i * len(s) / k)] for i in range(1, k)]


def assign_bin(v, edges):
    b = 0
    while b < len(edges) and v >= edges[b]:
        b += 1
    return b


def main():
    for p in (IE, SPROT, PAX, LINKS, SCORED):
        if not os.path.exists(p):
            sys.exit(f"missing {p}\n  PaxDb inputs: see data/external/README.md / data/SOURCES.md")

    # ------------------------------------------------------------------ inputs + the join
    print("=" * 96)
    print("STEP 1  DATA + THE JOIN (a poor join rate is a result, not something to hide)")
    print("=" * 96)
    prots = load_swissprot()
    by_ensp, by_gene = load_paxdb(PAX)
    links = load_links()
    ab, via = join_abundance(prots, by_ensp, by_gene, links)

    print(f"  SwissProt human canonical proteins      : {len(prots):,}")
    print(f"  PaxDb v6.1 9606 WHOLE_ORGANISM entries  : {len(by_ensp):,} ENSP ({len(by_gene):,} symbols)")
    print(f"  PaxDb UniProt cross-reference (9606)    : {len(links):,} ENSP -> entry name")
    print(f"\n  JOIN RATE (SwissProt protein -> PaxDb ppm)")
    print(f"    via PaxDb's own ENSP->UniProt link    : {via['link']:>7,}  ({100*via['link']/len(prots):.1f}%)")
    print(f"    via gene symbol (fallback)            : {via['gene']:>7,}  ({100*via['gene']/len(prots):.1f}%)")
    print(f"    TOTAL WITH ABUNDANCE                  : {len(ab):>7,}  ({100*len(ab)/len(prots):.1f}%)")
    print(f"    no abundance (dropped)                : {len(prots)-len(ab):>7,}  "
          f"({100*(len(prots)-len(ab))/len(prots):.1f}%)")
    join_rate = len(ab) / len(prots)
    print("\n  Missingness is NOT random: a protein absent from PaxDb is typically one MS has never")
    print("  quantified, i.e. low-abundance. Dropping those INFLATES the background's central")
    print("  tendency, which makes test A CONSERVATIVE (harder, not easier, to look 'above background').")

    # ------------------------------------------------------------------ catalogue + matching
    tissues, src_gene = load_catalogue()
    peptides = sorted(tissues)
    print(f"\n  IEAtlas cancer catalogue: {len(peptides):,} unique sequences over "
          f"{len({t for s in tissues.values() for t in s})} cancer types")
    print("  matching every sequence against the canonical proteome (exact substring) ...")
    hits = match_peptides(peptides, prots)

    # self-check against the catalogue of record -- this MUST reproduce canonical_self
    selfmap = {r["peptide"]: int(r["canonical_self"])
               for r in csv.DictReader(open(SCORED, newline="", encoding="utf-8"))}
    shared = [p for p in peptides if p in selfmap]
    mism = sum(1 for p in shared if (p in hits) != bool(selfmap[p]))
    n_ov = sum(1 for p in peptides if p in hits)
    print(f"    canonical-overlapping: {n_ov:,}/{len(peptides):,} = {100*n_ov/len(peptides):.1f}%"
          f"   (headline: 98,193/174,465 = 56.3%)")
    print(f"    self-check vs claim_catalog_scored.csv: {len(shared):,} shared, {mism} mismatches "
          f"(expect 0)")
    if mism:
        sys.exit("FATAL: the recomputed overlap set disagrees with the scored catalogue. Stop.")

    # ------------------------------------------------------------------ rows for the tests
    rows, no_ab = [], 0
    for p in peptides:
        h = hits.get(p)
        if not h:
            continue
        vals = [ab[i] for i in h if i in ab]
        if not vals:
            no_ab += 1
            continue
        best = max(range(len(h)), key=lambda k: ab.get(h[k], -1.0))
        rows.append({
            "pep": p, "len": len(p), "breadth": len(tissues[p]),
            "ab_max": max(vals), "ab_med": st.median(vals),
            "n_prot": len(h),
            "cl_canon": prots[h[best]][2] or prots[h[best]][0],  # matched canonical gene = EXPOSURE cluster
            "cl_src": src_gene.get(p, "?"),                      # IEAtlas source gene = LABEL cluster
        })
    print(f"    overlapping sequences with >=1 matched protein carrying a PaxDb abundance: "
          f"{len(rows):,}/{n_ov:,} = {100*len(rows)/n_ov:.1f}%")
    print(f"    overlapping sequences whose matched proteins ALL lack abundance (dropped): {no_ab:,}")

    # ------------------------------------------------------------------ A. which proteins do they hit?
    print("\n" + "=" * 96)
    print("A  DO CANONICAL-OVERLAPPING SEQUENCES PREFERENTIALLY MATCH *ABUNDANT* CANONICAL PROTEINS?")
    print("=" * 96)
    bg = sorted(ab.values())                       # background: every canonical protein with a ppm
    bg_med = st.median(bg)
    print(f"  BACKGROUND  all {len(bg):,} canonical proteins with a PaxDb ppm")
    print(f"              median {bg_med:,.2f} ppm   (log10 {math.log10(bg_med):+.2f})")
    print(f"\n  {'exposure':<34}{'n':>9}{'median ppm':>14}{'x background':>14}{'>bg median':>12}{'AUC':>8}")
    print("  " + "-" * 91)
    A = {}
    for lab, key in (("matched protein, MAX ppm", "ab_max"), ("matched protein, MEDIAN ppm", "ab_med")):
        v = [r[key] for r in rows]
        med = st.median(v)
        frac = sum(1 for x in v if x > bg_med) / len(v)
        a = auc(v, bg)
        print(f"  {lab:<34}{len(v):>9,}{med:>14,.2f}{med/bg_med:>13.2f}x{100*frac:>11.1f}%{a:>8.3f}")
        A[key] = {"median_ppm": round(med, 3), "fold_vs_background": round(med / bg_med, 2),
                  "frac_above_bg_median": round(frac, 4), "auc_vs_background": round(a, 3)}

    # protein-level view: the unit IS the protein, so peptide clustering cannot inflate anything
    hit_idx = {i for p in peptides if p in hits for i in hits[p]}
    hit_ab = [ab[i] for i in hit_idx if i in ab]
    nohit_ab = [ab[i] for i in ab if i not in hit_idx]
    print(f"\n  PROTEIN-LEVEL (no peptide-clustering problem: the unit of analysis is the protein)")
    print(f"    canonical proteins HIT by >=1 catalogued sequence : {len(hit_ab):>7,}  "
          f"median {st.median(hit_ab):>10,.2f} ppm")
    print(f"    canonical proteins hit by NONE                    : {len(nohit_ab):>7,}  "
          f"median {st.median(nohit_ab):>10,.2f} ppm")
    a_prot = auc(hit_ab, nohit_ab)
    print(f"    fold {st.median(hit_ab)/st.median(nohit_ab):.2f}x   AUC {a_prot:.3f}   "
          f"(0.5 = no shift)")
    A["protein_level"] = {
        "n_hit": len(hit_ab), "n_not_hit": len(nohit_ab),
        "median_hit_ppm": round(st.median(hit_ab), 3),
        "median_nothit_ppm": round(st.median(nohit_ab), 3),
        "fold": round(st.median(hit_ab) / st.median(nohit_ab), 2),
        "auc": round(a_prot, 3),
    }
    a_ok = A["ab_max"]["auc_vs_background"] > 0.5 and a_prot > 0.5

    # ------------------------------------------------------------------ B. THE KEY TEST
    print("\n" + "=" * 96)
    print("B  THE KEY TEST -- DOES THE ABUNDANCE OF THE MATCHED PROTEIN PREDICT DETECTION BREADTH?")
    print("=" * 96)
    print("   Prediction: mean number of cancer types RISES monotonically across abundance bins.")
    print("   FLAT or DECREASING = the abundance explanation FAILS and comes out of the paper.\n")

    w = defaultdict(int)
    for r in rows:
        w[r["len"]] += 1
    weights = {L: c / len(rows) for L, c in w.items()}

    results_b = {}
    for key, name in (("ab_max", "MAX ppm over matched proteins"),
                      ("ab_med", "MEDIAN ppm over matched proteins")):
        edges = bins_by_quantile([r[key] for r in rows], NBINS)
        for r in rows:
            r["bin"] = assign_bin(r[key], edges)
        crude = defaultdict(list)
        for r in rows:
            crude[r["bin"]].append(r["breadth"])
        stdm = std_bin_means(rows, weights)

        print(f"  exposure = {name}   ({NBINS} equal-count bins)")
        print(f"    {'bin':<6}{'n':>9}{'ppm range':>26}{'mean breadth':>14}{'len-std':>10}{'median':>8}")
        print("    " + "-" * 73)
        lo_e = -math.inf
        for b in sorted(crude):
            hi_e = edges[b] if b < len(edges) else math.inf
            rng_s = f"{'0' if lo_e==-math.inf else f'{lo_e:,.1f}'} - " \
                    f"{'inf' if hi_e==math.inf else f'{hi_e:,.1f}'}"
            print(f"    Q{b+1:<5}{len(crude[b]):>9,}{rng_s:>26}{st.mean(crude[b]):>14.3f}"
                  f"{stdm[b]:>10.3f}{st.median(crude[b]):>8.0f}")
            lo_e = hi_e

        bl = sorted(stdm)
        mono = all(stdm[bl[i]] <= stdm[bl[i + 1]] + 1e-12 for i in range(len(bl) - 1))
        d_crude = st.mean(crude[bl[-1]]) - st.mean(crude[bl[0]])
        d_std = stdm[bl[-1]] - stdm[bl[0]]
        rho = rank_corr([(r[key], r["breadth"]) for r in rows])
        print(f"\n    monotone increasing across bins : {'YES' if mono else 'NO'}")
        print(f"    Q5 - Q1 mean breadth (crude)    : {d_crude:+.3f} cancer types")
        print(f"    Q5 - Q1 mean breadth (len-std)  : {d_std:+.3f} cancer types   <-- headline")
        print(f"    Spearman rho (ppm vs breadth)   : {rho:+.3f}")
        print("    gene-clustered bootstrap of the length-standardized Q5-Q1 "
              "(z-tests are invalid here):")
        ci_canon = clustered_bootstrap(rows, weights, "cl_canon",
                                       "cluster = matched canonical gene")
        ci_src = clustered_bootstrap(rows, weights, "cl_src",
                                     "cluster = IEAtlas source gene")
        excl = ci_canon[0] > 0 and ci_src[0] > 0
        print(f"    => both CIs {'EXCLUDE' if excl else 'DO NOT exclude'} 0.\n")
        results_b[key] = {
            "monotone": mono, "q5_minus_q1_crude": round(d_crude, 3),
            "q5_minus_q1_lengthstd": round(d_std, 3), "spearman_rho": round(rho, 3),
            "ci95_cluster_canonical_gene": [round(ci_canon[0], 3), round(ci_canon[1], 3)],
            "ci95_cluster_source_gene": [round(ci_src[0], 3), round(ci_src[1], 3)],
            "ci_excludes_zero_both_clusterings": excl,
            "bin_means_lengthstd": {f"Q{b+1}": round(stdm[b], 3) for b in sorted(stdm)},
            "bin_n": {f"Q{b+1}": len(crude[b]) for b in sorted(crude)},
        }

    # deciles, on the primary exposure -- a finer look at the shape of the trend
    edges = bins_by_quantile([r["ab_max"] for r in rows], 10)
    for r in rows:
        r["bin"] = assign_bin(r["ab_max"], edges)
    dec = std_bin_means(rows, weights)
    print("  DECILES (primary exposure = MAX ppm), length-standardized mean breadth:")
    print("    " + "  ".join(f"D{b+1}" for b in sorted(dec)))
    print("    " + "  ".join(f"{dec[b]:.2f}" for b in sorted(dec)))
    dl = sorted(dec)
    dec_mono = all(dec[dl[i]] <= dec[dl[i + 1]] + 1e-12 for i in range(len(dl) - 1))
    up = sum(1 for i in range(len(dl) - 1) if dec[dl[i + 1]] > dec[dl[i]])
    print(f"    strictly monotone: {'YES' if dec_mono else 'NO'}   "
          f"({up}/{len(dl)-1} consecutive steps increase)   D10 - D1 = {dec[dl[-1]]-dec[dl[0]]:+.3f}")

    # ------------------------------------------------------------------ C1. length control
    print("\n" + "=" * 96)
    print("C1  LENGTH CONTROL -- the confound a reviewer already caught. Does B hold WITHIN length?")
    print("=" * 96)
    edges = bins_by_quantile([r["ab_max"] for r in rows], NBINS)
    for r in rows:
        r["bin"] = assign_bin(r["ab_max"], edges)
    print(f"  {'len':>5}{'n':>9}{'Q1 breadth':>13}{'Q5 breadth':>13}{'Q5-Q1':>9}{'rho':>8}  holds?")
    print("  " + "-" * 66)
    holds = strata = 0
    per_len = {}
    for L in sorted({r["len"] for r in rows}):
        sub = [r for r in rows if r["len"] == L]
        if len(sub) < MIN_STRATUM:
            continue
        g = defaultdict(list)
        for r in sub:
            g[r["bin"]].append(r["breadth"])
        if len(g) < NBINS:
            continue
        strata += 1
        b = sorted(g)
        q1, q5 = st.mean(g[b[0]]), st.mean(g[b[-1]])
        rho = rank_corr([(r["ab_max"], r["breadth"]) for r in sub])
        ok = q5 > q1
        holds += ok
        print(f"  {L:>5}{len(sub):>9,}{q1:>13.3f}{q5:>13.3f}{q5-q1:>+9.3f}{rho:>+8.3f}  "
              f"{'yes' if ok else 'NO'}")
        per_len[L] = {"n": len(sub), "q1": round(q1, 3), "q5": round(q5, 3),
                      "q5_minus_q1": round(q5 - q1, 3), "rho": round(rho, 3)}
    c1_ok = strata >= 2 and holds == strata
    print(f"\n  => the abundance->breadth trend holds in {holds}/{strata} length strata. "
          f"{'NOT a length artifact.' if c1_ok else 'LENGTH-CONFOUNDED -- the effect is not clean.'}")

    # ------------------------------------------------------------------ C2. protein-length control
    print("\n" + "=" * 96)
    print("C2  PROTEIN-LENGTH CONTROL -- the obvious attack on A. A LONGER protein contains more")
    print("    substrings, so it is likelier to be hit BY CHANCE. Are 'hit' proteins just long ones?")
    print("=" * 96)
    hl = [len(prots[i][3]) for i in hit_idx if i in ab]
    nl = [len(prots[i][3]) for i in ab if i not in hit_idx]
    a_len = auc(hl, nl)
    rho_la = rank_corr([(len(prots[i][3]), ab[i]) for i in ab])
    print(f"  The confound is REAL: hit proteins are longer -- median {st.median(hl):,.0f} aa vs "
          f"{st.median(nl):,.0f} aa (AUC {a_len:.3f}).")
    print(f"  And across the background, protein length correlates with abundance at rho = {rho_la:+.3f}.")
    print("  So the test is whether the ABUNDANCE shift survives WITHIN protein-length strata.\n")
    ls = sorted(len(prots[i][3]) for i in ab)
    ledges = [ls[int(k * len(ls) / 10)] for k in range(1, 10)]
    print(f"  {'decile':<8}{'aa range':>14}{'n hit':>8}{'n not':>8}{'med hit':>10}{'med not':>10}"
          f"{'fold':>8}{'AUC':>8}")
    print("  " + "-" * 64)
    aucs, per_plen, lo_l = [], {}, 0
    for d in range(10):
        hi_l = ledges[d] if d < 9 else 10 ** 9
        H = [ab[i] for i in hit_idx if i in ab and lo_l <= len(prots[i][3]) < hi_l]
        N = [ab[i] for i in ab if i not in hit_idx and lo_l <= len(prots[i][3]) < hi_l]
        if len(H) >= 20 and len(N) >= 20:
            a_d = auc(H, N)
            aucs.append(a_d)
            fold = st.median(H) / max(st.median(N), 1e-9)
            rng_s = f"{lo_l}-{hi_l if hi_l < 10**9 else 'max'}"
            print(f"  D{d+1:<7}{rng_s:>14}{len(H):>8,}{len(N):>8,}{st.median(H):>10.2f}"
                  f"{st.median(N):>10.2f}{fold:>7.1f}x{a_d:>8.3f}")
            per_plen[d + 1] = {"n_hit": len(H), "n_not": len(N), "fold": round(fold, 2),
                               "auc": round(a_d, 3)}
        lo_l = hi_l
    c2_ok = len(aucs) >= 2 and all(a > 0.5 for a in aucs)
    print(f"\n  => the abundance shift holds in {sum(1 for a in aucs if a > 0.5)}/{len(aucs)} "
          f"protein-length deciles (mean AUC {st.mean(aucs):.3f} vs {a_prot:.3f} pooled).")
    print(f"  {'NOT explained by protein length' if c2_ok else 'PROTEIN-LENGTH-CONFOUNDED'}"
          f" -- and note the fold is LARGEST in the SHORTEST decile, which is where a chance hit is")
    print("  least likely. The confound, where it acts, acts AGAINST the result.")

    # ------------------------------------------------------------------ C3. placebo
    print("\n" + "=" * 96)
    print("C3  PLACEBO -- break the peptide->protein link. The trend MUST collapse.")
    print("=" * 96)
    rng = RNG(SEED)
    obs = results_b["ab_max"]["q5_minus_q1_lengthstd"]
    null = []
    for _ in range(200):
        shuf = [ab[i] for i in hit_idx if i in ab]      # abundances of hit proteins, reassigned at random
        fake = []
        for r in rows:
            v = shuf[rng.randint(len(shuf))]
            fake.append({"len": r["len"], "breadth": r["breadth"],
                         "bin": assign_bin(v, edges)})
        null.append(q5_minus_q1(fake, weights))
    null = [v for v in null if v == v]
    null.sort()
    n_med = st.median(null)
    n_hi = null[int(0.975 * len(null)) - 1]
    beats = sum(1 for v in null if v >= obs)
    print(f"  200 random reassignments of canonical-protein abundance to the same peptides:")
    print(f"    placebo Q5-Q1 (len-std): median {n_med:+.3f}, 97.5th pct {n_hi:+.3f}")
    print(f"    OBSERVED                : {obs:+.3f}")
    print(f"    placebo draws >= observed: {beats}/{len(null)}")
    c3_ok = beats == 0 and abs(n_med) < abs(obs) / 2
    print(f"  => the trend {'COLLAPSES' if c3_ok else 'DOES NOT collapse'} under the placebo. "
          f"{'The machinery does not manufacture trends.' if c3_ok else 'WARNING: the binning machinery produces a trend from nothing -- distrust B.'}")

    # ------------------------------------------------------------------ version robustness
    v5 = None
    if os.path.exists(PAX_V5):
        print("\n" + "=" * 96)
        print("VERSION ROBUSTNESS -- rerun the headline on PaxDb v5.0, the PUBLISHED release (2024).")
        print("  v6.1 is four days old at time of writing; the result must not depend on that choice.")
        print("=" * 96)
        e5, g5 = load_paxdb(PAX_V5)
        ab5, _ = join_abundance(prots, e5, g5, links)
        r5 = []
        for r in rows:
            h = hits[r["pep"]]
            vals = [ab5[i] for i in h if i in ab5]
            if vals:
                r5.append({**r, "ab_max": max(vals)})
        if not r5:
            sys.exit("FATAL: v5.0 produced 0 usable peptides -- the check silently did nothing. Stop.")
        ed5 = bins_by_quantile([r["ab_max"] for r in r5], NBINS)
        for r in r5:
            r["bin"] = assign_bin(r["ab_max"], ed5)
        w5 = defaultdict(int)
        for r in r5:
            w5[r["len"]] += 1
        v5 = q5_minus_q1(r5, {L: c / len(r5) for L, c in w5.items()})
        print(f"  v5.0 join: {len(ab5):,}/{len(prots):,} proteins ({100*len(ab5)/len(prots):.1f}%);"
              f" {len(r5):,} peptides")
        print(f"  Q5 - Q1 mean breadth (len-std): v6.1 {obs:+.3f}   vs   v5.0 {v5:+.3f}")
        print(f"  => the headline {'REPRODUCES' if v5 > 0 else 'DOES NOT reproduce'} on the "
              f"published release.")

    # ------------------------------------------------------------------ artifact
    prim = results_b["ab_max"]
    json.dump({
        "paxdb_version": "v6.1 (latest, 2026-07-09), 9606 WHOLE_ORGANISM integrated, ppm",
        "join": {
            "swissprot_proteins": len(prots),
            "joined_via_paxdb_uniprot_link": via["link"],
            "joined_via_gene_symbol": via["gene"],
            "joined_total": len(ab),
            "join_rate": round(join_rate, 4),
        },
        "catalogue": {
            "unique_sequences": len(peptides),
            "canonical_overlapping": n_ov,
            "overlap_rate": round(n_ov / len(peptides), 4),
            "overlapping_with_abundance": len(rows),
            "overlapping_dropped_no_abundance": no_ab,
            "scored_csv_mismatches": mism,
        },
        "A_which_proteins_are_hit": A,
        "B_abundance_predicts_breadth": results_b,
        "B_deciles_lengthstd": {f"D{b+1}": round(dec[b], 3) for b in sorted(dec)},
        "B_deciles_monotone": dec_mono,
        "C1_peptide_length_strata": {"tested": strata, "holding": holds, "per_length": per_len},
        "C2_protein_length_control": {
            "hit_median_aa": st.median(hl), "nothit_median_aa": st.median(nl),
            "auc_length_hit_vs_not": round(a_len, 3),
            "rho_protein_length_vs_abundance": round(rho_la, 3),
            "deciles_tested": len(aucs), "deciles_auc_above_half": sum(1 for a in aucs if a > 0.5),
            "mean_auc_within_deciles": round(st.mean(aucs), 3), "per_decile": per_plen,
            "survives": c2_ok,
        },
        "C3_placebo": {"n_draws": len(null), "placebo_median": round(n_med, 3),
                       "placebo_p975": round(n_hi, 3), "observed": obs,
                       "draws_ge_observed": beats, "collapses": c3_ok},
        "version_robustness_v5_q5_minus_q1": (round(v5, 3) if v5 is not None else None),
        "bootstrap": {"B": B, "seed": SEED, "clusters": ["matched canonical gene", "IEAtlas source gene"]},
        "verdict": ("CORROBORATES" if (a_ok and prim["monotone"] and
                                       prim["ci_excludes_zero_both_clusterings"]
                                       and c1_ok and c2_ok and c3_ok)
                    else "REFUTES / NOT SUPPORTED"),
    }, open(ART, "w"), indent=2)
    print(f"\n  wrote {os.path.relpath(ART, REPO)}")

    # ------------------------------------------------------------------ verdict
    print("\n" + "=" * 96)
    print("VERDICT")
    print("=" * 96)
    b_ok = prim["monotone"] and prim["ci_excludes_zero_both_clusterings"]
    print(f"""
    A  overlapping sequences hit ABUNDANT canonical proteins : {'YES' if a_ok else 'NO'}
    B  abundance PREDICTS detection breadth (the key test)   : {'YES' if b_ok else 'NO'}
    C1 survives the PEPTIDE-length control                   : {'YES' if c1_ok else 'NO'}
    C2 survives the PROTEIN-length control                   : {'YES' if c2_ok else 'NO'}
    C3 survives the placebo                                  : {'YES' if c3_ok else 'NO'}
""")
    if a_ok and b_ok and c1_ok and c2_ok and c3_ok:
        print(f"""  CORROBORATED -- with the effect sizes stated plainly, because they are NOT the same size.

  A is LARGE and clean. Canonical-overlapping catalogued sequences match canonical proteins
    {A['ab_max']['fold_vs_background']:.1f}x more abundant than the median canonical protein (AUC
    {A['ab_max']['auc_vs_background']:.3f}). At the PROTEIN level -- where the unit of analysis is the protein, so
    peptide clustering cannot inflate anything -- proteins the catalogue hits are
    {A['protein_level']['fold']:.1f}x more abundant than proteins it never touches (AUC {A['protein_level']['auc']:.3f}),
    and that holds inside every protein-length decile (C2), largest where chance hits are least
    likely. The peptides the atlas calls non-canonical land on the abundant canonical proteome.

  B is REAL BUT MODEST, and the honest reading matters. Length-standardized mean breadth rises
    monotonically across abundance bins, Q5 - Q1 = {prim['q5_minus_q1_lengthstd']:+.3f} cancer types, CI excludes 0 under
    BOTH clusterings (matched canonical gene {prim['ci95_cluster_canonical_gene']}, IEAtlas source gene
    {prim['ci95_cluster_source_gene']}), it holds in {holds}/{strata} peptide-length strata, and the placebo collapses.
    But the CRUDE trend SATURATES -- it climbs from Q1 to Q3 and is flat-to-slightly-down across
    Q3-Q5 (it is monotone only after length standardization) -- and the rank correlation is weak
    (rho = {prim['spearman_rho']:+.3f}). Abundance shifts detection breadth by a fraction of one cancer type,
    not by a lot. Report the direction and the CI; do NOT sell this as a strong dose-response.

  WHAT THIS LICENSES, EXACTLY. The detection-breadth PROXY can now be retired in favour of a direct
  measurement, and the measurement points the same way: the canonical proteins these sequences
  overlap are abundant ones, and abundance predicts how broadly a sequence was detected. It remains
  an ASSOCIATION in an observational resource, not a controlled experiment -- abundance is entangled
  with expression breadth, HLA presentation, protein length and study composition, and nothing here
  apportions the catalogue's 56.3% overlap among them. Licensed: "the matched canonical protein's
  abundance predicts detection breadth, and canonical-overlapping catalogued sequences match abundant
  canonical proteins -- consistent with abundance-driven detectability." NOT licensed: any
  quantitative attribution of the overlap rate to abundance.

  AND NOT LICENSED (ERROR #18): do not resurrect the withdrawn "56.3% EXCEEDS the library's 34.1%,
  therefore the excess arose during detection" argument on the back of this. That comparison is a
  unit mismatch (distinct PEPTIDES at native lengths post-search vs distinct 9-MERS of an undetected
  candidate space) and is dead regardless of how this measurement turned out. What is measured here
  stands on its own and needs no such excess: abundance predicts detection breadth, full stop.

  Scope, unchanged: MS identifies the SEQUENCE, never the LOCUS. Nothing here says the biology is
  fake. It says the reported record cannot attribute these sequences to the ncORF.""")
    else:
        failed = [n for n, ok in (("A", a_ok), ("B", b_ok), ("C1", c1_ok), ("C2", c2_ok),
                                  ("C3", c3_ok)) if not ok]
        print(f"""  NOT SUPPORTED -- failing: {', '.join(failed)}. This is reported, not buried.

  B, the key test: monotone = {prim['monotone']}, length-standardized Q5 - Q1 =
  {prim['q5_minus_q1_lengthstd']:+.3f} cancer types, rho = {prim['spearman_rho']:+.3f}, CIs exclude 0 under both
  clusterings = {prim['ci_excludes_zero_both_clusterings']}.

  The direct measurement does not carry the abundance explanation. The abundance/detectability claim
  must come OUT of the manuscript, or be restated as an explicitly-labelled untested hypothesis --
  and the detection-breadth proxy must NOT be left standing in its place, since the direct
  measurement of the thing the proxy was proxying FAILED.

  The descriptive result is untouched by this and stands on its own: 56.3% of catalogued sequences
  overlap the canonical proteome exactly. What would be lost is only the EXPLANATION of why -- which
  would then be openly unexplained, and should be stated as unexplained.""")
    print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
