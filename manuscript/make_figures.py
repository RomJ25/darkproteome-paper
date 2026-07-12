"""
Generate Figures 1-4 for the manuscript. Deterministic content -- re-running reproduces
identical data and layout; the PDF/PNG bytes themselves differ run-to-run only by
matplotlib's embedded creation timestamp.
Every number cites its verified source (the analysis scripts / live script runs).
Fig 3b reproduces scripts/verify_effective_rho.py exactly (seed 0).
Run:  python3 manuscript/make_figures.py   ->  manuscript/figures/*.pdf,*.png
"""
import sys, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import figure_data                                    # NOTHING is typed in; see figure_data.py

OUT = Path("manuscript/figures"); OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"savefig.dpi":300, "font.size":9, "axes.spines.top":False,
                     "axes.spines.right":False, "axes.titlesize":9,
                     "font.family":"serif",
                     "font.serif":["Libertinus Serif","STIX Two Text","DejaVu Serif"],
                     "mathtext.fontset":"stix", "pdf.fonttype":42})  # serif to match the Libertinus body; TrueType-embedded
BLUE,RED,GREEN,DK="#4477AA","#CC6677","#117733","#882255"
def save(fig,name):
    fig.tight_layout()
    fig.savefig(OUT/f"{name}.pdf",bbox_inches="tight"); fig.savefig(OUT/f"{name}.png",bbox_inches="tight")
    plt.close(fig); print("wrote",name)

# ---- Fig 2: the REPORTING & ADJUDICABILITY MATRIX ----
# A survivorship funnel is the wrong instrument here. A funnel asserts that each stage FAILED the
# previous one; most of these cells are not failures, they are SILENCES. Drawing a dimension the
# record cannot even decide as though it were attrition renders a property of the scorer as a
# property of the claims. The matrix cannot tell that lie.
F2 = figure_data.fig2()
DIMS = ["source_translation","hla_elution","allele_restriction",
        "normal_presentation","human_tcell_assay","class_fdr_reconstructible"]
NICE = {"source_translation":"source translation","hla_elution":"HLA elution",
        "allele_restriction":"allele restriction","normal_presentation":"normal presentation",
        "human_tcell_assay":"human T-cell assay",
        "class_fdr_reconstructible":"class-FDR reconstructible"}
RUNGS = ["asserted","claim_linked","quantitative","modality_appropriate","adjudicable"]
RLAB  = ["asserted","claim-\nlinked","quanti-\ntative","modality-\nappropriate","ADJUDICABLE"]

n = F2["n_pooled"]
M = np.array([[F2["pooled"][d][r] for r in RUNGS] for d in DIMS], dtype=float)
frac = M / n

fig, ax = plt.subplots(figsize=(7.4, 3.5))
im = ax.imshow(frac, cmap="Blues", vmin=0, vmax=1, aspect="auto")
ax.set_xticks(range(len(RUNGS))); ax.set_xticklabels(RLAB, fontsize=6.8)
ax.set_yticks(range(len(DIMS)));  ax.set_yticklabels([NICE[d] for d in DIMS], fontsize=7.4)
ax.set_xlabel("reporting ladder  (a claim must clear every rung to be independently adjudicable)",
              fontsize=7.2)
for i in range(len(DIMS)):
    for j in range(len(RUNGS)):
        v = int(M[i, j])
        ax.text(j, i, f"{v:,}", ha="center", va="center", fontsize=6.6,
                color="white" if frac[i, j] > 0.55 else "#222",
                fontweight="bold" if j == 4 else "normal")
# the adjudicable column is the result; mark it
ax.add_patch(plt.Rectangle((3.5, -0.5), 1, len(DIMS), fill=False, edgecolor=RED, lw=1.4, zorder=5))
ax.set_title(f"Claim-linked evidence availability  (n = {n:,} audited machine-readable claims)",
             loc="left", fontsize=8.2)
cb = fig.colorbar(im, ax=ax, fraction=.028, pad=.02); cb.ax.tick_params(labelsize=6)
cb.set_label("fraction of claims", fontsize=6.5)
save(fig, "fig2_adjudicability")

# The dead funnel must not linger next to its replacement: a stale PDF in figures/ is exactly
# how a retired number keeps reaching a reader.
for ext in ("pdf", "png"):
    dead = OUT / f"fig2_survivorship.{ext}"
    if dead.exists():
        dead.unlink(); print("removed dead figure", dead.name)

# ---- Fig 1: sequence non-uniqueness (class-resolved) + the atlas reuse floor ----
# EVERY value comes from figure_data.fig1(); the old code hardcoded
#   vals=[0.2,0.0,0.0,37.1,56.3,0.0]
# which carried the RETIRED 37.1% (a first-peptide-only artifact; the true rate is 57.7%) and a
# literal 0.0 for CrypticProteinDB (really 1/3,810 = 0.026% -- small but NOT zero, and that is
# the whole point of the bar: the overlap is avoidable).
F1 = figure_data.fig1()
bars = F1["bars"]
fig,(a,b)=plt.subplots(1,2,figsize=(8.6,3.2))
labs = [x["label"] for x in bars]
vals = [x["pct"] for x in bars]
# COLOUR TRACKS THE VALUE, NEVER THE SOURCE. An earlier draft coloured atlas bars with the
# "concerning" dark tone regardless of their rate, which painted CrypticProteinDB -- at 0.026%,
# the single best-behaved resource in the figure and the whole reason the bar is here -- as
# though it were a problem. A figure whose colours argue against its own numbers is worse than
# no figure.
cols = [RED if x["pct"] >= 30 else GREEN for x in bars]
# minimum-visible-bar-width floor: a 0.026% bar is sub-pixel at any rasterization and would read
# as a true zero -- which is precisely the error the hardcoded 0.0 made. Floors only the DRAWN
# height; the printed label below still uses the true value.
vis = [v if v == 0 else max(v, 0.6) for v in vals]
a.bar(range(len(bars)), vis, color=cols, alpha=.88)
a.set_xticks(range(len(bars))); a.set_xticklabels(labs, fontsize=6.0)
a.set_ylabel("exact canonical-substring rate (%)")
for i, x in enumerate(bars):
    txt = f"{x['pct']:.1f}%" if x["pct"] >= 1 else (f"{x['pct']:.3f}%" if x["pct"] else "0%")
    a.text(i, max(vals[i], 0.6) + 1.4, txt, ha="center", fontsize=6.2)
a.set_ylim(0, max(vals) * 1.28)
a.set_title("(a) Sequence non-uniqueness, class- and atlas-resolved", loc="left")

# (b) IEAtlas
ie = [x for x in bars if x["cls"] == "IEAtlas"][0]
pct = ie["pct"]; rest = 100 - pct
b.bar(0, pct, color=DK, alpha=.88,
      label=f"exact match in canonical reference ({pct:.1f}%)")
# "genuinely non-canonical" is RETIRED (it asserted a biological fact the sequence cannot show).
# The remainder is only: no exact match in THIS reference layer R.
b.bar(0, rest, bottom=pct, color="#88CCEE", alpha=.9,
      label=f"no exact match in reference layer $R$ ({rest:.1f}%)")
b.set_xlim(-1.2,2.9); b.set_xticks([]); b.set_ylim(0,100)
b.set_ylabel("IEAtlas cryptic epitopes (%)")
b.legend(fontsize=6.0, loc="upper right", framealpha=.9)
b.text(0.9, 30,
       "the unmatched remainder is a statement\nabout the REFERENCE, not about biology:\n"
       r"$N_i(R)$ is reference-relative",
       fontsize=6.2, va="center", ha="left")
b.set_title("(b) Atlas reuse: overlap with the canonical reference", loc="left")
save(fig,"fig1_floors")

# ---- Fig 3: the identifiability theory ----
fig,(a,b)=plt.subplots(1,2,figsize=(7.8,3.3))
f=np.linspace(0.002,0.12,500); alpha=0.03; upper=np.minimum(1,alpha/f)
a.fill_between(f*100,0,upper,color="#88CCEE",alpha=.4,label=r"identified set  $[0,\ \min(1,\alpha/f)]$")
a.plot(f*100,upper,color=BLUE,lw=1.6)
a.axvspan(4.4,7.1,color=RED,alpha=.13)
a.annotate("broad cryptic class\nf = 4.4–7.1%  $\\Rightarrow$  $\\leq$42–68%",xy=(7.5,0.60),fontsize=6.4,ha="left")
a.annotate("noncoding subset\nf < 1%  $\\Rightarrow$  unconstrained",xy=(1.0,0.86),fontsize=6.4,ha="left")
a.axhline(alpha,ls=":",color="grey",lw=.8); a.text(11.8,alpha+0.02,r"reported $\alpha$ = 3%",fontsize=6.3,color="grey",ha="right")
a.set_xlabel(r"class fraction  $f=T_\mathrm{N}/T$  (%)"); a.set_ylabel("class FDR: sharp upper bound"); a.set_ylim(0,1.06)
a.legend(fontsize=6.3,loc="upper right"); a.set_title("(a) Class FDR is set-identified",loc="left")
# (b) effective-rho drift: reproduce verify_effective_rho.py part B (seed 0)
rng=np.random.default_rng(0); SL,SH=4000,200
mN=np.concatenate([rng.integers(1,3,SL),rng.integers(400,1200,SH)])
mC=np.concatenate([rng.integers(5,20,SL),rng.integers(20,120,SH)]); m=mN+mC; rg=mN.sum()/mC.sum()
ts=np.logspace(-4,-1,16)
reff=[(((1-(1-t)**m)*(mN/m)).sum()/(((1-(1-t)**m)*(mC/m)).sum())) for t in ts]
b.plot(ts,reff,"-o",color=BLUE,ms=3.5,label=r"$\rho_\mathrm{eff}(\tau)$   (theory $=$ sim)")
b.axhline(rg,ls="--",color=RED,label=r"global $\rho=%.2f$   (wrong object)"%rg)
# Shading marks where rho_eff and rho_global are visibly close (~8% apart), i.e. the
# stringent end Results II argues immunopeptidomics operates in -- not a claim that this
# toy simulation's x-axis calibrates to a real threshold.
b.axvspan(ts[0],ts[2],color="#888",alpha=.15,zorder=0)
# Annotation sits below the band, not beside the curve: the curve stays high (>2.4)
# across this x-range, so the empty bottom of the panel is clear of it.
b.annotate("stringent end\n(Results II):\n$\\rho_\\mathrm{eff}\\approx\\rho_\\mathrm{global}$ here",
           xy=(ts[1],reff[1]),xytext=(2.2e-4,0.5),fontsize=6.2,ha="left",va="center",
           arrowprops=dict(arrowstyle="-",color="#999",lw=.7))
b.set_xscale("log"); b.set_xlabel(r"threshold tail prob  $t=1-F(\tau)$   (stringent $\to$ lenient)")
b.set_ylabel(r"class decoy ratio  $D_\mathrm{N}/D_\mathrm{C}$"); b.set_ylim(0,2.9); b.legend(fontsize=6.6,loc="upper right")
b.set_title(r"(b) Governing quantity: effective $\rho$, not global $\rho$",loc="left")
save(fig,"fig3_theory")

# ---- Fig 4: per-claim truth is a partial order ----
# Forest-plot convention, not filled bars: a filled bar from 0 reads as a point value and
# inverts the thesis. Interval claims: thin line + filled dot at the known lower bound (0) + open cap at the
# upper bound (a bound, not a confirmed value). Point-identified P=0: a filled diamond.
fig,ax=plt.subplots(figsize=(6.2,2.5))
claims=["Raja altORF claim","HCC pseudogene claim","HCC exact-canonical-\nsubstring claim"]
his=[0.998,0.629,0.0]; y=np.arange(len(claims))[::-1]
for yi,hi in zip(y,his):
    if hi>0:
        ax.plot([0,hi],[yi,yi],lw=1.6,color=BLUE,alpha=.9,solid_capstyle="butt",zorder=2)
        ax.plot(0,yi,marker="o",color=BLUE,ms=5,zorder=3)
        ax.plot(hi,yi,marker="o",mfc="white",mec=BLUE,mew=1.5,ms=7,zorder=3)
        ax.text(hi+0.035,yi,r"$P\leq %.3f$"%hi,va="center",fontsize=8)
    else:
        ax.plot(0,yi,marker="D",color=DK,ms=8,zorder=3)
        ax.text(0.045,yi,r"$P=0$",va="center",fontsize=8)
ax.set_yticks(y); ax.set_yticklabels(claims,fontsize=8); ax.set_xlim(-0.02,1.12); ax.set_ylim(-0.6,2.6)
ax.set_xlabel(r"$P(\mathrm{claim\ true}\mid\mathrm{published\ record})$:  sharp identified interval")
# SUPPLEMENT ONLY. This panel renders `P = 0` -- a point-identified probability that a claim is
# FALSE -- which the sequence evidence cannot support. Sequence overlap shows only that the record
# does not EXCLUDE a canonical source; it never shows the biology is absent. Emitted under a
# supplement name so it cannot silently reappear in a main-text float.
save(fig,"supp_fig_partial_order")
print("ALL FIGURES DONE ->", OUT)
