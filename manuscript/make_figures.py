"""
Generate Figures 1-4 for the manuscript. Reproducible (re-run -> identical figures).
Every number cites its verified source (the analysis scripts / live script runs).
Fig 3b reproduces scripts/verify_effective_rho.py exactly (seed 0).
Run:  python3 manuscript/make_figures.py   ->  manuscript/figures/*.pdf,*.png
"""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

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

# ---- Fig 2: re-verifiability survivorship (reference_model.py) ----
# Linear scale (not log): on a log axis the 54% first-step cut and the 8.5% second-step
# cut render as nearly the same bar-length change as the 940x allele-assignment cliff,
# anti-correlated with true attrition. Linear also lets the two terminal zero counts be
# plotted at x=0 instead of footnoted.
# Monochrome, not a red highlight on the last 3 stages: red is reserved for the
# "concerning finding" color elsewhere (Fig 1a/1b).
stages=["all unique non-canonical peptides","not canonical-self (provenance-clean)",
        "absent from normal (tumor-restricted)","allele-assigned (presentation)",
        "source-translation substantiated","human T-cell validated"]
counts=[180650,82453,75471,80,0,0]
fig,ax=plt.subplots(figsize=(6.6,3.3)); y=np.arange(len(stages))[::-1]
xmax=180650
# minimum-visible-bar-width floor: the raw 80-peptide (0.044%) bar is sub-pixel at any
# rasterization resolution and renders visually identical to the true-zero rows below it.
# Floors only the DRAWN bar length for nonzero counts; every label below still uses the
# true `counts`.
vis_counts=[c if c==0 else max(c,xmax*0.004) for c in counts]
ax.barh(y,vis_counts,color=BLUE,alpha=.88)
ax.set_yticks(y); ax.set_yticklabels(stages,fontsize=8)
def pstr(c):
    p=100*c/180650
    if c==0: return "0 (0%)"
    if abs(p-100)<1e-6: return f"{c:,} (100%)"
    return f"{c:,} ({p:.1f}%)" if p>=1 else f"{c:,} ({p:.2g}%)"
for yi,c in zip(y,counts):
    ax.text(c+xmax*0.017,yi,pstr(c),va="center",fontsize=7.5)
# marginal-independence note: the funnel's sequential order isn't load-bearing for the
# 0% headline — translation alone is reported for ≈0% of the corpus (Results I).
ax.text(xmax*0.10,1.0,"translation is independently $\\approx$0% of the\ncorpus — not solely a downstream effect\nof allele-assignment (Results I)",
        fontsize=6.3,color="#444",va="center",ha="left")
ax.set_xlim(0,xmax*1.30); ax.set_xlabel("peptides")  # no in-figure title — the LaTeX caption describes it
save(fig,"fig2_survivorship")

# ---- Fig 1: non-novelty (class-resolved) + IEAtlas reuse contamination ----
fig,(a,b)=plt.subplots(1,2,figsize=(8.2,3.2))
labs=["altORF\n(Raja)","lncRNA\n(Raja)","pseudo.\n(Raja)","pseudo.\n(HCC)","IEAtlas\n(atlas)","Cryptic-\nProteinDB"]
vals=[0.2,0.0,0.0,37.1,56.3,0.0]; cols=[GREEN,GREEN,GREEN,RED,DK,GREEN]
# same minimum-visible-bar-width floor as Fig 2 above: the 0.2% altORF bar is thin enough
# that its 1px antialiased edge blends with the axis spine instead of the intended fill
# color. Floors only the drawn height; the label loop below still uses the true `vals`.
vis_vals=[v if v==0 else max(v,0.5) for v in vals]
a.bar(range(6),vis_vals,color=cols,alpha=.88)
a.set_xticks(range(6)); a.set_xticklabels(labs,fontsize=6.6); a.set_ylabel("exact canonical-substring rate (%)")
for i,v in enumerate(vals): a.text(i,v+1.2,f"{v:g}%",ha="center",fontsize=7)
a.set_title("(a) Non-novelty floor — class- and atlas-resolved",loc="left")
# (b) IEAtlas: 174,465 cancer cryptic epitopes
# Leader-line callout, not a stacked sub-bar, for the 9.1% fact: on this 0-100%-of-whole-bar
# axis a sub-bar sized to "9.1% of the 43.7% remainder" would render only a few points
# tall, mismatching its own label.
b.bar(0,56.3,color=DK,alpha=.88,label="canonical-self by sequence (56.3%)")
b.bar(0,43.7,bottom=56.3,color="#88CCEE",alpha=.9,label="genuinely non-canonical (43.7%)")
b.set_xlim(-1.2,2.6); b.set_xticks([]); b.set_ylim(0,100); b.set_ylabel("IEAtlas cryptic epitopes (%)")
b.legend(fontsize=6.3,loc="upper right",framealpha=.9)  # clear of the bar (bar spans x=[-0.4,0.4]); callouts sit below it
b.annotate("9.1% of this 43.7% slice\nappears on normal tissue",
           xy=(0.42,70),xytext=(0.85,58),fontsize=6.4,va="center",ha="left",
           arrowprops=dict(arrowstyle="-",color="#666",lw=.7))
# No leader line into the bar: 16/43 is the HCC cohort's own canonical-self pseudogene
# count (panel a), not an IEAtlas statistic -- naming the population directly in the text
# avoids misattributing it to this bar.
b.text(0.85,25,"pseudogene specificity (HCC cohort,\npanel a): 16/43 directly on normal\n(HLA Ligand Atlas)",
       fontsize=6.4,va="center",ha="left")
b.set_title("(b) IEAtlas reuse: canonical-self + a normal floor",loc="left")
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
a.set_xlabel(r"class fraction  $f=T_\mathrm{N}/T$  (%)"); a.set_ylabel("class FDR — sharp upper bound"); a.set_ylim(0,1.06)
a.legend(fontsize=6.3,loc="upper right"); a.set_title("(a) Class FDR is set-identified",loc="left")
# (b) effective-rho drift — reproduce verify_effective_rho.py part B (seed 0)
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
ax.set_xlabel(r"$P(\mathrm{claim\ true}\mid\mathrm{published\ record})$  —  sharp identified interval")
save(fig,"fig4_partial_order")
print("ALL FIGURES DONE ->", OUT)
