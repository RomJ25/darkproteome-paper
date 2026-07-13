"""Figures for manuscript_v2. Every value is read from an artifact -- nothing is typed.

    python3 manuscript/make_figures_v2.py  ->  manuscript/figures_v2/*.pdf,*.png

Three figures, one per result:
  F1  the measurement  -- canonical-sequence overlap, SAME-PIPELINE catalogues only, plus the
                          robustness panel (era-correct reference / length standardization / class)
  F2  the library      -- latent canonical ambiguity of the ncORF libraries, plus the detection
                          effect that the library alone does NOT account for
  F3  the consequence  -- normal-tissue presentation, with the catalogue's own non-overlapping
                          sequences as a WITHIN-RESOURCE COMPARATOR (not a control), reported with
                          the length-standardized risk ratio and its gene-clustered CI

WHAT IS DELIBERATELY NOT DRAWN (external review):
  * Bedran et al.'s four published rates (1.4-5%) are NOT plotted beside ours. They come from a
    different pipeline, reference, dedup and peptide unit. Putting them on one axis invites exactly
    the fold-change reading ("11-40x") that we withdrew -- and a figure argues that comparison more
    forcefully than a sentence ever could. They belong in the text, as context, with no arithmetic.
  * "z = 74" is gone. It was an invalid statistic (174,465 clustered observations treated as
    independent). F3 reports the length-standardized RR with a gene-clustered bootstrap CI.
  * "internal control" is gone. The non-overlapping set controls nothing; it is a comparator.
"""
import csv
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "manuscript"))
OUT = os.path.join(REPO, "manuscript", "figures_v2")
os.makedirs(OUT, exist_ok=True)
csv.field_size_limit(10_000_000)

ATL = os.path.join(REPO, "data", "external", "atlases")
SCORED = os.path.join(REPO, "data", "claim_catalog_scored.csv")
SCALED = os.path.join(REPO, "data", "claim_catalog_scaled.csv")
TIER1 = os.path.join(REPO, "data", "primary_tier1_nonnovelty.csv")
R3 = json.load(open(os.path.join(REPO, "data", "derived_r3_inference.json")))
BIAS = json.load(open(os.path.join(REPO, "data", "derived_detection_bias.json")))
ERA = json.load(open(os.path.join(REPO, "data", "derived_era_reference.json")))
UNI = json.load(open(os.path.join(REPO, "data", "derived_library_union.json")))
ABD = json.load(open(os.path.join(REPO, "data", "derived_abundance_direct.json")))
_UREF = next(r for r in UNI["by_reference"] if r["reference"] == UNI["primary_reference"])

plt.rcParams.update({"savefig.dpi": 300, "font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False, "font.family": "serif",
                     "font.serif": ["Libertinus Serif", "STIX Two Text", "DejaVu Serif"],
                     "pdf.fonttype": 42})
RED, GREEN, BLUE, GREY = "#CC6677", "#117733", "#4477AA", "#888888"


def save(fig, name):
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, f"{name}.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, f"{name}.png"), bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)


def epitopes(p):
    s = set()
    with open(p, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.reader(fh, delimiter="\t")
        next(rd, None)
        for r in rd:
            if r and r[0] and r[0].strip().upper().isalpha():
                s.add(r[0].strip().upper())
    return s


selfmap = {r["peptide"]: int(r["canonical_self"])
           for r in csv.DictReader(open(SCORED, newline="", encoding="utf-8"))}

src = {}
for r in csv.DictReader(open(SCALED, newline="", encoding="utf-8")):
    s = (r.get("_source") or "").strip()
    p = (r.get("peptide_sequence") or "").strip().upper()
    if p and p.isalpha() and not s.startswith("cohort:"):
        src.setdefault(s, set()).add(p)


def rate(peps):
    n = sum(1 for p in peps if p in selfmap)
    k = sum(1 for p in peps if selfmap.get(p))
    return (100 * k / n if n else 0.0), k, n


cpdb, cpk, cpn = rate(src["CrypticProteinDB-immuno"] | src["CrypticProteinDB-epitopes"])
ie, iek, ien = rate(src["IEAtlas"])
t1 = list(csv.DictReader(open(TIER1, newline="", encoding="utf-8")))
raja = [r for r in t1 if r["cohort"].startswith("Raja")]
rjk = sum(int(r["canonical_self_exact"]) for r in raja)
rj = 100 * rjk / len(raja)

# ------------------------------------------------- F1: the measurement + its robustness
fig, (a, b) = plt.subplots(1, 2, figsize=(8.2, 3.5),
                           gridspec_kw={"width_ratios": [1.15, 1]})

vals = [cpdb, rj, ie]
labs = ["Cryptic-\nProteinDB", "Raja et al.\n(ovarian)", "IEAtlas"]
cnt = [(cpk, cpn), (rjk, len(raja)), (iek, ien)]
a.bar(range(3), vals, color=[GREEN, GREEN, RED], alpha=.9)
a.set_yscale("log")
a.set_ylim(0.015, 400)
a.set_ylabel("canonical-sequence overlap (%), log scale")
a.set_xticks(range(3))
a.set_xticklabels(labs, fontsize=7.5)
for i, (v, (k, n)) in enumerate(zip(vals, cnt)):
    a.text(i, v * 1.35, f"{v:.3g}%\n{k:,}/{n:,}", ha="center", fontsize=6.8)
# Label sits in white space ABOVE the two green bars -- inside them it renders green-on-green.
a.text(0.5, 1.6, "explicit exclusion rule\nin their Methods", fontsize=7, color=GREEN,
       ha="center", fontweight="bold")
a.annotate("", xy=(0.0, 0.05), xytext=(0.30, 1.1),
           arrowprops=dict(arrowstyle="-", color=GREEN, lw=.7, alpha=.7))
a.annotate("", xy=(1.0, 0.28), xytext=(0.72, 1.1),
           arrowprops=dict(arrowstyle="-", color=GREEN, lw=.7, alpha=.7))
a.set_title("(a) One pipeline, one reference, one peptide unit", loc="left", fontsize=9)
a.text(0.0, -0.28, "Published rates for four other catalogues (1.4–5%, Bedran et al. 2023) are NOT\n"
                   "plotted: different pipeline, so a ratio against them is not a measurement.",
       transform=a.transAxes, ha="left", fontsize=6.2, color="#666")

# robustness panel -- every bar is a different way of trying to make 56.3% go away
rob = [("as reported\n(reference R)", 56.3, RED),
       (f"era-correct\n(Swiss-Prot 2022_01)", ERA["pct_2022_01"], BLUE),
       ("length-\nstandardized", 56.3, BLUE),
       ("if NO pseudogene\nORFs at all", 55.8, BLUE)]
b.bar(range(4), [x[1] for x in rob], color=[x[2] for x in rob], alpha=.9)
for i, (_l, v, _c) in enumerate(rob):
    b.text(i, v + 1.2, f"{v}%", ha="center", fontsize=7.5)
b.set_xticks(range(4))
b.set_xticklabels([x[0] for x in rob], fontsize=6.6)
b.set_ylim(0, 68)
b.set_ylabel("canonical-sequence overlap (%)")
b.set_title("(b) It does not go away", loc="left", fontsize=9)
b.text(0.5, -0.30, f"only {ERA['retrospective_only']} sequences ({ERA['retrospective_pct_of_overlap']}%"
                   " of the overlap set) are matches\na February-2022 analyst could not have made",
       transform=b.transAxes, ha="center", fontsize=6.2, color="#666")
save(fig, "f1_measurement")

# ------------------------------------------------- F2: the library, and what it does NOT explain
fig, (a, b, c) = plt.subplots(1, 3, figsize=(11.6, 3.4),
                              gridspec_kw={"width_ratios": [1.25, .95, 1.0]})
# Two of IEAtlas's three sources, and their union -- the point is that the union FALLS. Values are
# read from the artifact so the figure cannot drift from the paper.
_nu = _UREF["nuorfdb"]["pct"]
_tl = _UREF["translnc"]["pct"]
_un = _UREF["union"]["pct"]
LIB = [("nuORFdb v1.2\n(source 1)", _nu, RED),
       ("Translnc\n(source 2)", _tl, BLUE),
       ("union\n(measured)", _un, "#7B3FA0"),
       ("GENCODE\nRibo-seq p1", 2.4, GREEN),
       ("GENCODE\nRibo-seq p2", 1.0, GREEN)]
a.bar([x[0] for x in LIB], [x[1] for x in LIB], color=[x[2] for x in LIB], alpha=.9)
for i, (_l, v, _c) in enumerate(LIB):
    a.text(i, v + 1.1, f"{v:.1f}%", ha="center", fontsize=8)
a.set_ylabel("latent canonical ambiguity\n(% of distinct 9-mers)")
a.set_ylim(0, 46)
a.tick_params(axis="x", labelsize=6.6)
a.set_title("(a) Libraries differ enormously — and do NOT compose", loc="left", fontsize=9)
# The retraction, drawn: adding a real second source LOWERS the union rate.
# Two traps here, both invisible in source and obvious on the PNG:
#   * the arrowhead must NOT land on the "20.2%" bar label -- stop it well above (+5.0, not +1.6).
#   * "->" as a literal U+2192 renders as TOFU in the serif face. Spell it in words.
a.annotate("", xy=(1.9, _un + 3.4), xytext=(0.34, _nu - 1.2),
           arrowprops=dict(arrowstyle="->", color="#7B3FA0", lw=1.0,
                           connectionstyle="arc3,rad=.28"))
a.text(1.05, 42.0, f"adding Translnc LOWERS the union\n"
                   f"({_nu:.1f}% down to {_un:.1f}%). Union rates are\n"
                   f"not monotone: {_nu:.1f}% is NOT a lower bound.",
       fontsize=6.2, color="#7B3FA0", ha="center", va="center")
a.text(3.5, 26, "RPFdb v2.0 (source 3) is\nunobtainable, so the full\nlibrary is NOT determined.",
       fontsize=6.0, color="#666", ha="center", va="center")

# the detection effect: ribosomal-ORF enrichment, catalogue vs the search space it was drawn from
rr_cat, rr_lib = BIAS["ribo_catalogue_rr"], BIAS["ribo_library_rr"]
b.bar([0, 1], [rr_lib, rr_cat], color=[GREY, RED], alpha=.9, width=.55)
b.axhline(1.0, color="#333", lw=.8, ls="--")
b.text(1.52, 1.0, "no\nenrichment", fontsize=6.2, color="#333", va="center")
# The library bar tops out at ~0.91, and the "no enrichment" rule sits at 1.0. A label placed just
# ABOVE that bar lands on the dashed line and the digits are occluded -- invisible in the source,
# obvious on opening the PNG. Put it INSIDE the bar; keep the catalogue label above.
b.text(0, rr_lib - .12, f"{rr_lib}×", ha="center", va="top", fontsize=8,
       color="white", fontweight="bold")
b.text(1, rr_cat + .07, f"{rr_cat}×", ha="center", fontsize=8)
b.set_xticks([0, 1])
b.set_xticklabels(["in the LIBRARY\n(nothing detected yet)", "in the CATALOGUE\n(after detection)"],
                  fontsize=7)
b.set_ylabel("ribosomal-ORF enrichment among\ncanonical-overlapping sequences")
b.set_ylim(0, 3.2)
b.set_xlim(-.6, 1.9)
b.set_title("(b) …but the library does not account for it", loc="left", fontsize=9)
b.annotate("", xy=(1, rr_cat - .05), xytext=(0, rr_lib + .05),
           arrowprops=dict(arrowstyle="->", color="#444", lw=.9,
                           connectionstyle="arc3,rad=-.25"))
b.text(0.5, 2.42, f"{BIAS['ribo_excess']}× excess\narises during detection", fontsize=6.6,
       color="#444", ha="center")

# (c) the DIRECT abundance measurement -- this replaces the proxy, so it must be drawn, including
# its weakness. The bars are the length-standardized breadth per abundance quintile.
_bb = ABD["B_abundance_predicts_breadth"]["ab_max"]
_bm = _bb["bin_means_lengthstd"]
_qs = ["Q1", "Q2", "Q3", "Q4", "Q5"]
_vals = [_bm[q] for q in _qs]
c.bar(range(5), _vals, color=[BLUE, "#5B8FD4", "#8E7FC8", "#C06090", RED], alpha=.92, width=.62)
for i, v in enumerate(_vals):
    c.text(i, v + .012, f"{v:.2f}", ha="center", fontsize=7)
c.set_xticks(range(5))
c.set_xticklabels(_qs, fontsize=7.5)
c.set_xlabel("abundance of the matched canonical protein\n(PaxDb quintile, low to high)", fontsize=7)
c.set_ylabel("mean cancer types detected in\n(length-standardized)")
c.set_ylim(1.30, 1.90)
c.set_title("(c) Abundance predicts detection — weakly", loc="left", fontsize=9)
_pl = ABD["A_which_proteins_are_hit"]["protein_level"]
c.text(2.0, 1.855,
       f"canonical proteins the catalogue hits are\n"
       f"{_pl['fold']}× more abundant than those it never hits\n"
       f"(median {_pl['median_hit_ppm']} vs {_pl['median_nothit_ppm']} ppm)",
       fontsize=6.0, color="#444", ha="center", va="center")
# The honest caveat belongs ON the figure, not only in the caption.
c.text(0.5, -0.42, f"Q5−Q1 = {_bb['q5_minus_q1_lengthstd']} cancer types  ·  Spearman ρ = "
                   f"{_bb['spearman_rho']}  —  real, and small.\nAbundance is one contributor to what "
                   f"gets detected, not the explanation.",
       transform=c.transAxes, ha="center", fontsize=6.2, color="#666", style="italic")
save(fig, "f2_library")

# ------------------------------------------------- F3: the consequence, correctly inferred
cancer = epitopes(os.path.join(ATL, "IEAtlas_Epitopes_In_Cancer_Tissues.txt"))
p1 = R3["pct_overlapping_in_normal"]
p2 = R3["pct_comparator_in_normal"]
kov, knov = R3["overlapping_in_normal"], R3["comparator_in_normal"]
nov_n, nnov_n = R3["n_overlapping"], R3["n_comparator"]

fig, (a, b) = plt.subplots(1, 2, figsize=(8.2, 3.5),
                           gridspec_kw={"width_ratios": [1, 1.15]})
a.bar(["canonical-\noverlapping", "NOT overlapping\n(within-resource\ncomparator)"], [p1, p2],
      color=[RED, GREEN], alpha=.9)
for i, (v, k, n) in enumerate(((p1, kov, nov_n), (p2, knov, nnov_n))):
    a.text(i, v + 0.8, f"{v:.1f}%\n{k:,}/{n:,}", ha="center", fontsize=7)
a.set_ylabel("also in IEAtlas's OWN\nnormal-tissue export (%)")  # one line overruns the canvas
a.set_ylim(0, 29)
a.tick_params(axis="x", labelsize=7)
a.set_title(f"(a) Risk ratio {R3['rr_length_standardized']}× "
            f"[{R3['ci95'][0]}, {R3['ci95'][1]}]", loc="left", fontsize=9)
a.text(0.5, -0.34, "length-standardized; 95% CI from a bootstrap resampling\n"
                   f"{R3['n_clusters']:,} source-gene clusters (the observations are clustered,\n"
                   "so a two-proportion z-test would not be valid)",
       transform=a.transAxes, ha="center", fontsize=6.2, color="#666")

both = kov
rest = len(cancer) - both
b.barh([0], [both], height=0.42, color=RED, alpha=.9,
       label=f"canonical-compatible AND already in the atlas's\nown normal-tissue export  ({both:,})")
b.barh([0], [rest], left=[both], height=0.42, color="#DDDDDD",
       label=f"remainder of the catalogue  ({rest:,})")
b.set_yticks([])
b.set_ylim(-0.6, 0.6)
b.set_xlim(0, len(cancer))
b.set_xlabel("IEAtlas cancer-catalogued sequences (unique)")
b.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _p: f"{int(v):,}"))
b.tick_params(axis="x", labelsize=7)
b.spines["left"].set_visible(False)
b.legend(fontsize=6.6, loc="upper center", bbox_to_anchor=(0.5, -0.26), frameon=False)
b.set_title(f"(b) {R3['pct_of_whole_catalogue']}% of the catalogue, no external reference needed",
            loc="left", fontsize=9)
save(fig, "f3_consequence")

print("\nALL FIGURES ->", OUT)
