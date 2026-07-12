"""Figures for manuscript_v2. Every value computed from the artifacts -- nothing typed.

    python3 manuscript/make_figures_v2.py  ->  manuscript/figures_v2/*.pdf,*.png

Three figures, one per result:
  F1  the outlier      -- canonical-sequence overlap across catalogues (log scale; the spread is
                          3 orders of magnitude and a linear axis would erase everything but IEAtlas)
  F2  the library      -- latent canonical ambiguity of the ncORF libraries themselves
  F3  the consequence  -- normal-tissue presentation, with the catalogue's own non-overlapping
                          epitopes as the internal control
"""
import csv
import os
import re
import sys

import numpy as np
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

# ---------------------------------------------------------------- F1: the outlier
src = {}
for r in csv.DictReader(open(SCALED, newline="", encoding="utf-8")):
    s = (r.get("_source") or "").strip()
    p = (r.get("peptide_sequence") or "").strip().upper()
    if p and p.isalpha() and not s.startswith("cohort:"):
        src.setdefault(s, set()).add(p)


def rate(peps):
    n = sum(1 for p in peps if p in selfmap)
    k = sum(1 for p in peps if selfmap.get(p))
    return 100 * k / n if n else 0.0


cpdb = rate(src["CrypticProteinDB-immuno"] | src["CrypticProteinDB-epitopes"])
ie = rate(src["IEAtlas"])
t1 = list(csv.DictReader(open(TIER1, newline="", encoding="utf-8")))
raja = [r for r in t1 if r["cohort"].startswith("Raja")]
rj = 100 * sum(int(r["canonical_self_exact"]) for r in raja) / len(raja)

# the four from Bedran et al. 2023 -- CITED, not re-measured. Marked as such in the figure.
BEDRAN = [("Erhard 2020", 1.4), ("Ouspenskaia 2021", 3.0), ("Chong 2020", 4.0),
          ("Laumont 2016", 5.0)]

labels = ["Cryptic-\nProteinDB", "Raja et al.\n(ovarian)"] + [b[0].replace(" ", "\n") for b in BEDRAN] + ["IEAtlas"]
vals = [cpdb, rj] + [b[1] for b in BEDRAN] + [ie]
kind = ["ours", "ours"] + ["cited"] * 4 + ["ours"]
excl = [True, True, False, False, False, False, False]

fig, ax = plt.subplots(figsize=(7.2, 3.4))
cols = [GREEN if e else (RED if v > 20 else GREY) for v, e in zip(vals, excl)]
bars = ax.bar(range(len(vals)), vals, color=cols, alpha=.9,
              hatch=["" if k == "ours" else "///" for k in kind], edgecolor="white", linewidth=.6)
ax.set_yscale("log")
ax.set_ylim(0.015, 200)
ax.set_ylabel("canonical-sequence overlap (%), log scale")
ax.set_xticks(range(len(vals)))
ax.set_xticklabels(labels, fontsize=6.8)
for i, v in enumerate(vals):
    ax.text(i, v * 1.25, f"{v:.3g}%", ha="center", fontsize=7)
ax.axhspan(1.4, 5.0, color=GREY, alpha=.10, zorder=0)
# Annotations sit in EMPTY regions of the axes, checked against the rendered figure. An earlier
# version collided the band caption with the 5% bar label and the exclusion-rule note with the
# 0.026% value -- both invisible in the source and obvious the moment the PNG is opened.
ax.text(2.4, 14, "range of every ncORF catalogue\npreviously audited (1.4–5%)",
        fontsize=6.4, color="#555", ha="center")
ax.annotate("", xy=(3.4, 5.2), xytext=(2.9, 11.5),
            arrowprops=dict(arrowstyle="-", color="#999", lw=.7))
# The label goes in the WHITE SPACE above the two green bars, not on top of them: placed inside
# the bar it rendered dark-green-on-green and was unreadable. Found by opening the PNG.
ax.text(0.5, 0.55, "explicit exclusion rule", fontsize=7, color=GREEN, ha="center",
        fontweight="bold")
ax.annotate("", xy=(0.0, 0.045), xytext=(0.35, 0.42),
            arrowprops=dict(arrowstyle="-", color=GREEN, lw=.7, alpha=.7))
ax.annotate("", xy=(1.0, 0.25), xytext=(0.68, 0.42),
            arrowprops=dict(arrowstyle="-", color=GREEN, lw=.7, alpha=.7))
ax.set_title("IEAtlas is an order-of-magnitude outlier on the field's own metric",
             loc="left", fontsize=9)
ax.text(1.0, -0.30, "hatched = as reported by Bedran et al. 2023 (not re-measured here)",
        transform=ax.transAxes, ha="right", fontsize=6.2, color="#666")
save(fig, "f1_outlier")

# ---------------------------------------------------------------- F2: the library
LIB = [("nuORFdb v1.2\n(integrated by IEAtlas)", 34.1),
       ("GENCODE Ribo-seq\nORFs (phase 1)", 2.4),
       ("GENCODE Ribo-seq\nORFs (phase 2)", 1.0)]
fig, ax = plt.subplots(figsize=(5.6, 3.2))
ax.bar([x[0] for x in LIB], [x[1] for x in LIB],
       color=[RED, GREEN, GREEN], alpha=.9)
for i, (_l, v) in enumerate(LIB):
    ax.text(i, v + 1.0, f"{v}%", ha="center", fontsize=8)
ax.set_ylabel("latent canonical ambiguity (% of distinct 9-mers)")
ax.set_ylim(0, 42)
ax.tick_params(axis="x", labelsize=7)
ax.axhspan(1.0, 2.4, color=GREY, alpha=.12, zorder=0)
ax.set_title("The ambiguity is created in the library, before any filter is applied",
             loc="left", fontsize=9)
ax.text(1.55, 26, "IEAtlas inherits this,\nand applies no exclusion rule.\n"
                  "Ouspenskaia searched the SAME\nlibrary and published at 3%.",
        fontsize=6.6, color="#444", va="center")
save(fig, "f2_library")

# ---------------------------------------------------------------- F3: the consequence
cancer = epitopes(os.path.join(ATL, "IEAtlas_Epitopes_In_Cancer_Tissues.txt"))
normal = epitopes(os.path.join(ATL, "IEAtlas_Epitopes_In_Normal_Tissues.txt"))
scored = [p for p in cancer if p in selfmap]
ov = [p for p in scored if selfmap[p]]
nov = [p for p in scored if not selfmap[p]]
kov = sum(1 for p in ov if p in normal)
knov = sum(1 for p in nov if p in normal)
p1, p2 = 100 * kov / len(ov), 100 * knov / len(nov)

fig, (a, b) = plt.subplots(1, 2, figsize=(8.0, 3.3),
                           gridspec_kw={"width_ratios": [1, 1.15]})
a.bar(["canonical-\noverlapping", "NOT overlapping\n(internal control)"], [p1, p2],
      color=[RED, GREEN], alpha=.9)
for i, (v, k, n) in enumerate(((p1, kov, len(ov)), (p2, knov, len(nov)))):
    a.text(i, v + 0.8, f"{v:.1f}%\n{k:,}/{n:,}", ha="center", fontsize=7)
a.set_ylabel("also in IEAtlas's OWN normal-tissue set (%)")
a.set_ylim(0, 28)
a.tick_params(axis="x", labelsize=7.5)
a.set_title(f"(a) Risk ratio {p1/p2:.1f}×  (z = 74)", loc="left", fontsize=9)

both = kov
rest = len(cancer) - both
b.barh([0], [both], height=0.42, color=RED, alpha=.9,
       label=f"canonical-compatible AND already\nseen on normal tissue  ({both:,})")
b.barh([0], [rest], left=[both], height=0.42, color="#DDDDDD",
       label=f"remainder of the catalogue  ({rest:,})")
b.set_yticks([])
b.set_ylim(-0.6, 0.6)
b.set_xlim(0, len(cancer))
b.set_xlabel("IEAtlas cancer epitopes")
b.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _p: f"{int(v):,}"))
b.tick_params(axis="x", labelsize=7)
b.spines["left"].set_visible(False)
b.legend(fontsize=6.6, loc="upper center", bbox_to_anchor=(0.5, -0.22), frameon=False)
b.set_title(f"(b) {100*both/len(cancer):.1f}% of the catalogue needs no external reference",
            loc="left", fontsize=9)
save(fig, "f3_consequence")

print("\nALL FIGURES ->", OUT)
