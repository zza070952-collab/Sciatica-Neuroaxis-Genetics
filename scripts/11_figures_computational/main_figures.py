from __future__ import annotations

import math
import os
import textwrap
from pathlib import Path
import warnings

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle, Circle
from scipy import stats

warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = Path(os.environ.get("SCIATICA_PROJECT_ROOT", Path.cwd()))
INPUT = Path(os.environ.get("PHASE8_SOURCE_DATA", ROOT / "data/computational_source_data"))
PHASE = Path(os.environ.get("PHASE8_OUT", ROOT / "results/phase8_JTM_core_final"))
DRG_RECURRENCE = Path(os.environ.get(
    "PHASE8_DRG_RECURRENCE",
    PHASE / "03_drg_recurrence_fix/DRG_expression_positive_recurrence.tsv",
))
OUT = PHASE / "06_main_figures"
SD = PHASE / "10_source_data/main_figure_tsv"
PREVIEW = PHASE / "14_validation/main_figure_preview"
AUDIT = PHASE / "14_validation/main_figure_text_audit"
for p in (OUT, SD, PREVIEW, AUDIT):
    p.mkdir(parents=True, exist_ok=True)

SCI = "#1F6F78"; LDH = "#B9852B"; TX = "#B65346"; MAPK = "#3E78A8"
FGFR = "#D3A332"; PDPR = "#766B88"; GFPT = "#8B7E91"; NAVY = "#19364D"
GREY = "#697680"; LIGHT = "#E7EDF0"; PALE = "#F7F9FA"; DARK = "#263740"
GENE_COL = {"TXNL1": TX, "MAPK3": MAPK, "FGFR3": FGFR, "PDPR": PDPR, "GFPT1": GFPT}
TISSUE_COL = {
    "Nerve_Tibial": "#0072B2", "Brain_Spinal_cord_cervical_c-1": "#CC79A7",
    "Whole_Blood": "#D55E00", "Muscle_Skeletal": "#009E73",
    "Cells_Cultured_fibroblasts": "#E69F00",
}
TISSUE_LABEL = {
    "Nerve_Tibial": "Tibial nerve", "Brain_Spinal_cord_cervical_c-1": "Spinal cord",
    "Whole_Blood": "Whole blood", "Muscle_Skeletal": "Skeletal muscle",
    "Cells_Cultured_fibroblasts": "Fibroblast",
}
CELL_COL = {
    "Excitatory neuron": "#4C78A8", "Inhibitory neuron": "#72B7B2", "Astrocyte": "#F2CF5B",
    "Oligodendrocyte": "#B279A2", "OPC": "#FF9DA6", "Microglia": "#E45756",
    "Endothelial": "#59A14F", "Pericyte": "#9D755D", "Ependymal": "#BAB0AC",
}
CELL_SHORT = {
    "Ambiguous broad identity": "Ambiguous", "Immune_non-microglia": "Immune (other)",
    "Meningeal/fibroblast-like": "Meningeal", "Motor/cholinergic neuron": "Motor/cholinergic",
    "Excitatory neuron": "Excitatory", "Inhibitory neuron": "Inhibitory",
    "Oligodendrocyte": "Oligodendro.",
}


def setup():
    mpl.rcParams.update({
        "font.family": "Liberation Sans", "font.size": 7.2, "axes.titlesize": 8.0,
        "axes.labelsize": 7.3, "xtick.labelsize": 6.7, "ytick.labelsize": 6.7,
        "axes.linewidth": 0.65, "axes.spines.top": False, "axes.spines.right": False,
        "legend.fontsize": 6.4, "pdf.fonttype": 42, "ps.fonttype": 42,
        "svg.fonttype": "none", "savefig.facecolor": "white", "figure.facecolor": "white",
    })


_CACHE: dict[tuple[int, str], pd.DataFrame] = {}


SHEET_MAP = {
    1: {"A_H_numeric_foundation": "A_H_numeric_foundation"},
    2: {"A_GWAS": "A_GWAS", "B_loci": "B_loci", "B_QQ": "C_QQ", "C_D_LDSC": "D_E_LDSC", "E_rg": "F_rg"},
    3: {"A": "A", "B": "B", "C": "C", "D": "D", "E": "E"},
    4: {"A": "A_FUSION", "B": "C_top_multisnp", "C": "D_joint", "D": "E_TXNL1_locus", "E": "F_coloc", "G": "G_context_coverage", "B_heatmap": "B_candidate_heatmap"},
    5: {"A": "A_hierarchy", "F": "F_donor", "G": "G_meta", "H": "H_LODO", "I": "I_concordance", "J": "J_controls", "K": "K_TXNL1"},
    6: {"A": "A_UMAP", "B": "C_composition", "C": "E_spinal_expression", "E": "F_DRG_expression", "markers": "B_markers", "assoc": "D_celltype_association", "registry": "G_H_registry"},
}


def read_sheet(n: int, sheet: str) -> pd.DataFrame:
    key = (n, sheet)
    if key not in _CACHE:
        _CACHE[key] = pd.read_excel(INPUT / f"Source_Data_Figure{n}.xlsx", sheet_name=sheet)
    return _CACHE[key].copy()


def load_fig(n: int, name: str) -> pd.DataFrame:
    if n == 5 and name in "BCDE":
        maps = read_sheet(5, "B_E_maps")
        donors = ["Donor19", "Donor43", "Donor45", "Donor47"]
        return maps[maps["donor"].eq(donors["BCDE".index(name)])].copy()
    return read_sheet(n, SHEET_MAP[n][name])


def load_sup(n: int) -> pd.DataFrame:
    raise RuntimeError("Phase 7C figures use locked per-figure Source Data workbooks")


def panel(ax, letter, subtitle=None, x=-0.12, y=1.08, title_pad=8):
    text = ax.text(x, y, letter, transform=ax.transAxes, ha="left", va="bottom", fontsize=10.2,
                   fontweight="bold", color=NAVY, clip_on=False)
    text.set_gid(f"panel_letter_{letter}")
    if subtitle:
        title = ax.set_title(subtitle, loc="left", fontsize=8.0, fontweight="semibold", color=DARK, pad=title_pad)
        title.set_gid(f"panel_title_{letter}")


def tidy(ax, grid="y"):
    if grid:
        ax.grid(axis=grid, color="#DCE3E7", lw=.45, zorder=0)
    ax.tick_params(length=2.2, width=.6, pad=1.5)


def shrink_figure_elements(fig, text_factor=.86, marker_factor=.80, min_text=5.35):
    """Apply a conservative, layout-only scale reduction to dense figures."""
    for txt in fig.findobj(mpl.text.Text):
        if not txt.get_visible() or not txt.get_text().strip():
            continue
        role = _text_role(txt)
        if role == "panel_letter":
            txt.set_fontsize(max(8.8, txt.get_fontsize() * .90))
        elif role == "panel_title":
            txt.set_fontsize(max(6.3, txt.get_fontsize() * text_factor))
        else:
            txt.set_fontsize(max(min_text, txt.get_fontsize() * text_factor))
    for collection in fig.findobj(mpl.collections.PathCollection):
        sizes = collection.get_sizes()
        if len(sizes):
            collection.set_sizes(np.maximum(2.0, sizes * marker_factor))
    for line in fig.findobj(mpl.lines.Line2D):
        if line.get_marker() not in (None, "", "None"):
            line.set_markersize(max(2.2, line.get_markersize() * math.sqrt(marker_factor)))


def strip(fig, text, y, color=SCI, h=.023):
    fig.add_artist(Rectangle((.018, y), .964, h, transform=fig.transFigure, fc=color, ec="none", zorder=20))
    fig.text(.028, y+h/2, text, color="white", va="center", ha="left", fontsize=7.3,
             fontweight="bold", zorder=21)


def _text_role(text):
    gid = text.get_gid() or ""
    if gid.startswith("panel_letter"):
        return "panel_letter"
    if gid.startswith("panel_title"):
        return "panel_title"
    ax = getattr(text, "axes", None)
    if ax is not None and text in getattr(ax, "texts", []):
        return "annotation"
    return "axis_or_legend"


def audit_matplotlib_text(fig, n):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    texts = [x for x in fig.findobj(mpl.text.Text) if x.get_visible() and x.get_text().strip()]
    items = []
    for idx, txt in enumerate(texts):
        try:
            box = txt.get_window_extent(renderer=renderer)
        except Exception:
            continue
        items.append((idx, txt, box))
    rows = []
    for i, (idx_a, a, ba) in enumerate(items):
        for idx_b, b, bb in items[i + 1:]:
            inter = mpl.transforms.Bbox.intersection(ba, bb)
            if inter is None or inter.width * inter.height <= 1.5:
                continue
            if a.axes is b.axes and _text_role(a) == _text_role(b) == "axis_or_legend" and a.get_text() == b.get_text():
                continue
            rows.append({
                "figure": f"Figure{n}", "object_a_id": idx_a, "object_b_id": idx_b,
                "role_a": _text_role(a), "role_b": _text_role(b),
                "text_a": a.get_text().replace("\n", " / "), "text_b": b.get_text().replace("\n", " / "),
                "intersection_px2": round(inter.width * inter.height, 2),
            })
    fw, fh = fig.canvas.get_width_height()
    for idx, txt, box in items:
        if box.x0 < -0.5 or box.y0 < -0.5 or box.x1 > fw + 0.5 or box.y1 > fh + 0.5:
            rows.append({
                "figure": f"Figure{n}", "object_a_id": idx, "object_b_id": "BOUNDARY",
                "role_a": _text_role(txt), "role_b": "figure_boundary",
                "text_a": txt.get_text().replace("\n", " / "), "text_b": "OUTSIDE",
                "intersection_px2": 0,
            })
    columns = ["figure","object_a_id","object_b_id","role_a","role_b","text_a","text_b","intersection_px2"]
    pd.DataFrame(rows, columns=columns).to_csv(AUDIT / f"Figure{n}_matplotlib_text_audit.tsv", sep="\t", index=False)


def save(fig, n, dpi=600):
    audit_matplotlib_text(fig, n)
    fig.savefig(OUT / f"Figure{n}_final.pdf")
    fig.savefig(OUT / f"Figure{n}_final.svg")
    fig.savefig(OUT / f"Figure{n}_final.tiff", dpi=dpi,
                pil_kwargs={"compression": "tiff_lzw"})
    for scale, preview_dpi in ((100, 300), (75, 225), (50, 150)):
        fig.savefig(PREVIEW / f"Figure{n}_final_{scale}pct.png", dpi=preview_dpi)
    plt.close(fig)


def write_source(n, name, d):
    p = SD / f"Figure{n}"
    p.mkdir(parents=True, exist_ok=True)
    d.to_csv(p / f"{name}.tsv", sep="\t", index=False)


def card(ax, x, y, w, h, title, lines, color, icon=None):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=.008,rounding_size=.015",
                                fc="white", ec=color, lw=1.0))
    ax.add_patch(Rectangle((x, y+h-.035), w, .035, fc=color, ec="none"))
    ax.text(x+.014, y+h-.0175, title, color="white", fontsize=6.9, fontweight="bold", va="center")
    yy = y+h-.055
    for left, right in lines:
        ax.text(x+.014, yy, left, color=GREY, fontsize=5.6, va="top")
        ax.text(x+w-.014, yy, right, color=DARK, fontsize=6.3, fontweight="bold", va="top", ha="right")
        yy -= .032
    if icon == "chromosome":
        ax.plot([x+.025,x+.045],[y+.025,y+.045],color=color,lw=4,solid_capstyle="round")
        ax.plot([x+.045,x+.025],[y+.025,y+.045],color=color,lw=4,solid_capstyle="round")


def cumulative_x(d, chr_col="chr", pos_col="position"):
    z = d.copy(); z[chr_col] = pd.to_numeric(z[chr_col], errors="coerce")
    z[pos_col] = pd.to_numeric(z[pos_col], errors="coerce"); z = z.dropna(subset=[chr_col,pos_col])
    z[chr_col] = z[chr_col].astype(int)
    maxs = z.groupby(chr_col)[pos_col].max().reindex(range(1,23)).fillna(0)
    offsets = maxs.shift(fill_value=0).cumsum()
    z["cum_x"] = z[pos_col] + z[chr_col].map(offsets)
    mids = offsets + maxs/2
    return z, mids, offsets


def fig1():
    f = plt.figure(figsize=(7.0866, 5.25))
    ax = f.add_axes([.02,.02,.96,.96]); ax.set(xlim=(0,1), ylim=(0,1)); ax.axis("off")
    strip(f, "PHENOTYPE AND GENETIC FOUNDATION", .926, SCI, .028)
    card(ax,.025,.685,.25,.205,"A  PHENOTYPE HIERARCHY",[("Strict sciatica","28,094 cases"),("Controls","347,768"),("LDH structural reference","36,573 cases")],SCI,"chromosome")
    card(ax,.295,.685,.30,.205,"B  GWAS FOUNDATION",[("Released variants","21.3 million"),("Common high-quality","8.28 million"),("HapMap3 variants","1.217 million"),("Observed h2","0.0382"),("Sciatica–LDH rg","0.9339")],SCI)
    # miniature Manhattan
    g = load_fig(2,"A_GWAS").dropna(subset=["x","neglog10p"])
    iax=f.add_axes([.625,.70,.34,.16]); iax.scatter(g.x,g.neglog10p,s=.3,c=np.where(g.chr.astype(int)%2,SCI,"#74A9B0"),rasterized=True)
    iax.axhline(-math.log10(5e-8),ls="--",lw=.55,c=TX); iax.set_xticks([]); iax.set_yticks([0,8,16]); iax.tick_params(labelsize=6.0,pad=1); iax.set_ylabel("-log10 P",fontsize=6.0,labelpad=2); iax.text(.02,.95,"Strict sciatica GWAS",transform=iax.transAxes,ha="left",va="top",fontsize=6.6,fontweight="bold",bbox=dict(fc="white",ec="none",alpha=.8,pad=.7)); tidy(iax,None)
    strip(f, "GENE AND REGULATORY REFINEMENT", .625, PDPR, .028)
    card(ax,.025,.395,.25,.195,"C  MAGMA + FUSION",[("MAGMA tested / FDR","18,341 / 640"),("FUSION models / FDR","40,437 / 398"),("Joint-retained models","11")],PDPR)
    # micro lollipop/TWAS
    iax=f.add_axes([.36,.41,.135,.145]); vals=[640,398,11]; iax.barh([2,1,0],vals,color=[SCI,"#43A39A",TX]); iax.set_yticks([2,1,0],["MAGMA","TWAS","Joint"],fontsize=5.8); iax.tick_params(axis="y",pad=3); iax.set_xticks([]); iax.set_xlim(0,800); iax.spines[:].set_visible(False)
    for yv,v in zip([2,1,0],vals): iax.text(v+14,yv,f"{v}",va="center",ha="left",fontsize=6.0,fontweight="bold",color=DARK)
    card(ax,.515,.395,.255,.195,"D  FULL-CIS COLOCALIZATION",[("Pre-spatial genes","20"),("Usable tests","118"),("Strong ABF contexts","15"),("SuSiE-supported genes","4")],PDPR)
    # cis signal glyph
    xs=np.linspace(.792,.91,45); y1=np.exp(-((xs-.846)/.018)**2); y2=.88*np.exp(-((xs-.849)/.021)**2)
    ax.plot(xs,.46+.08*y1,c=SCI,lw=1.2); ax.plot(xs,.43+.07*y2,c=TX,lw=1.2); ax.text(.853,.555,"aligned cis signals",fontsize=5.6,ha="center",color=GREY)
    ax.add_patch(Rectangle((.934,.385),.04,.22,fc="#FFF4F2",ec=TX,lw=1)); ax.text(.954,.495,"GENETICS\nFREEZE",rotation=90,ha="center",va="center",fontsize=6.0,fontweight="bold",color=TX)
    strip(f, "ADULT NEUROAXIS LOCALIZATION AND FUNCTIONAL VALIDATION", .335, "#406A72", .028)
    card(ax,.025,.075,.24,.22,"F  ADULT LUMBAR SPATIAL",[("GSE222322","20 sections"),("Biological replicates","4 donors"),("","section → donor → meta")],"#406A72")
    sp=load_fig(5,"B"); iax=f.add_axes([.06,.085,.15,.075]); q=iax.scatter(sp.spatial_x,sp.spatial_y,c=sp.z,s=1.5,cmap="RdBu_r",vmin=-4,vmax=4,rasterized=True); iax.axis("off")
    card(ax,.285,.075,.25,.22,"G  SINGLE-NUCLEUS CONTEXT",[("Spinal cord","68,175 / 9 donors"),("DRG","1,837 / 6 prep."),("Unique DRG donors","5")],"#406A72")
    u=load_fig(6,"A").sample(4000,random_state=7); iax=f.add_axes([.325,.085,.15,.075]); iax.scatter(u.UMAP1,u.UMAP2,s=.3,c=u.cell_type.map(CELL_COL).fillna("#B8C0C4"),rasterized=True); iax.axis("off")
    for x in [.27,.58,.765]: ax.add_patch(FancyArrowPatch((x,.79),(x+.015,.79),arrowstyle="-|>",mutation_scale=9,color=GREY,lw=.8))
    ax.add_patch(FancyArrowPatch((.536,.185),(.553,.185),arrowstyle="-|>",mutation_scale=8,color=TX,lw=1.0))
    write_source(1,"A_H_numeric_foundation",pd.DataFrame({"metric":["sciatica_cases","controls","LDH_cases","released_variants","common_HQ","HapMap3","observed_h2","rg","MAGMA_genes","MAGMA_FDR","FUSION_models","FUSION_FDR","joint_retained","pre_spatial_genes","coloc_tests","strong_ABF","SuSiE_genes","spatial_sections","spatial_donors","spinal_nuclei","spinal_donors","DRG_nuclei","DRG_preparations","DRG_unique_donors"],"value":[28094,347768,36573,21325018,8280937,1217311,.0382,.9339,18341,640,40437,398,11,20,118,15,4,20,4,68175,9,1837,6,5]}))
    save(f,1,600)


def fig2():
    f=plt.figure(figsize=(7.0866,5.55)); gs=f.add_gridspec(2,4,height_ratios=[1.75,1],hspace=.53,wspace=.68)
    A=load_fig(2,"A_GWAS"); A,mids,offs=cumulative_x(A,"chr","position")
    ax=f.add_subplot(gs[0,:3]); panel(ax,"A","Strict sciatica genome-wide association",x=-.055)
    for ch,g in A.groupby("chr"): ax.scatter(g.cum_x,g.neglog10p,s=.45,c=SCI if ch%2 else "#7CAEB2",alpha=.78,rasterized=True)
    ax.axhline(-math.log10(5e-8),c=TX,ls="--",lw=.75); ax.axvspan(offs.loc[6]+25e6,offs.loc[6]+34e6,color="#D9DDDF",alpha=.55)
    ax.set_xticks(mids.values); ax.set_xticklabels(range(1,23),fontsize=5.8); ax.set(ylabel="-log10(P)",xlabel="Chromosome"); ax.margins(x=.018); tidy(ax,None)
    for tick in ax.get_xticklabels()[-3:]: tick.set_fontsize(4.2)
    ax.text(.99,.95,"28,094 cases  |  347,768 controls",transform=ax.transAxes,ha="right",va="top",fontsize=6.7,color=GREY)
    # annotate descriptive candidate loci only
    cand={"FGFR3":(4,1808663),"MAPK3":(16,30000000),"TXNL1":(18,56600000)}
    for gene,(ch,pos) in cand.items():
        g=A[A.chr.eq(ch)]; r=g.iloc[(g.position-pos).abs().argsort()[:1]]
        if len(r):
            rr=r.iloc[0]; ax.annotate(gene,(rr.cum_x,rr.neglog10p),xytext=(0,9),textcoords="offset points",ha="center",fontsize=6,c=GENE_COL[gene],fontweight="bold",arrowprops=dict(arrowstyle="-",lw=.5,color=GENE_COL[gene]))
    loci=load_fig(2,"B_loci").copy()
    axb=f.add_subplot(gs[0,3]); panel(axb,"B","Descriptive lead\nclusters",x=-.18,title_pad=8)
    loci=loci.drop_duplicates("cluster_id").sort_values("lead_p").head(13); loci["nlp"]=-np.log10(pd.to_numeric(loci.lead_p))
    axb.hlines(np.arange(len(loci)),0,loci.nlp,color=LIGHT,lw=2); axb.scatter(loci.nlp,np.arange(len(loci)),c=SCI,s=16,zorder=3)
    cluster_short=[str(c).split("_")[-1] for c in loci.cluster_id]
    axb.set_yticks(np.arange(len(loci)),[f"chr{int(x)} · {c}" for x,c in zip(loci.chr,cluster_short)],fontsize=6.0); axb.tick_params(axis="y",pad=2); axb.invert_yaxis(); axb.set_xlabel("Lead -log10(P)"); tidy(axb,"x")
    q=load_fig(2,"B_QQ"); ax=f.add_subplot(gs[1,0]); panel(ax,"C","QQ calibration",x=-.16)
    ax.plot(q.expected_neglog10p,q.observed_neglog10p,c=SCI,lw=1); lim=max(q.expected_neglog10p.max(),q.observed_neglog10p.max()); ax.plot([0,lim],[0,lim],ls="--",c=GREY,lw=.7)
    n=8280937; ranks=np.linspace(1,n,250); lo=-np.log10(stats.beta.ppf(.975,ranks,n-ranks+1)); hi=-np.log10(stats.beta.ppf(.025,ranks,n-ranks+1)); exp=-np.log10(ranks/(n+1)); order=np.argsort(exp); ax.fill_between(exp[order],lo[order],hi[order],color=LIGHT,alpha=.65,lw=0)
    ax.set(xlabel="Expected -log10(P)",ylabel="Observed -log10(P)"); ax.text(.04,.93,"lambdaGC = 1.19\nLDSC intercept = 1.093",transform=ax.transAxes,va="top",fontsize=5.9); tidy(ax,None)
    ld=load_fig(2,"C_D_LDSC"); obs=ld[ld.scale.eq("observed")]
    ax=f.add_subplot(gs[1,1]); panel(ax,"D","Observed-scale h²",x=-.16)
    order=["SCIATICA","LDH"]; o=obs.set_index("trait").loc[order]
    for i,(trait,color) in enumerate(zip(order,[SCI,LDH])):
        ax.errorbar(i,o.loc[trait,"h2"],yerr=1.96*o.loc[trait,"se"],fmt="o",color=color,ecolor=color,capsize=3,ms=5)
    ax.set_xticks([0,1],["Sciatica","LDH"]); ax.set_ylabel("SNP h2 (95% CI)"); top=max(o.h2+1.96*o.se)+.014; ax.set_ylim(min(0.034, min(o.h2-1.96*o.se)-.004),top); tidy(ax)
    for i,r in enumerate(o.itertuples()): ax.text(i,r.h2+1.96*r.se+.0015,f"Z={r.z:.1f}",ha="center",va="bottom",fontsize=6.1)
    li=ld[ld.scale.eq("liability_sensitivity") & ld.trait.eq("SCIATICA")].sort_values("population_prevalence_assumption")
    ax=f.add_subplot(gs[1,2]); panel(ax,"E","Liability-scale sensitivity",x=-.16)
    ax.errorbar(li.population_prevalence_assumption,li.h2,yerr=1.96*li.se,fmt="o-",c=SCI,ms=3,capsize=2); ax.set(xlabel="Population prevalence K",ylabel="Liability h2"); tidy(ax)
    rg=load_fig(2,"E_rg").iloc[0]; ax=f.add_subplot(gs[1,3]); panel(ax,"F","Shared genetic\narchitecture",x=-.16,title_pad=8)
    ax.errorbar(rg.rg,0,xerr=1.96*rg.se,fmt="o",ms=6,c=SCI,capsize=3); ax.axvline(0,c=GREY,lw=.6); ax.set(xlim=(0,1.08),yticks=[],xlabel="Sciatica–LDH rg (95% CI)"); ax.text(.03,.82,f"rg={rg.rg:.4f}\n95% CI {rg.rg-1.96*rg.se:.3f}–{rg.rg+1.96*rg.se:.3f}",transform=ax.transAxes,fontsize=6.2,fontweight="bold"); tidy(ax,"x")
    for name,d in [("A_GWAS",A),("B_loci",loci),("C_QQ",q),("D_E_LDSC",ld),("F_rg",load_fig(2,"E_rg"))]: write_source(2,name,d)
    # The same-release/shared-control and structural-reference boundary is stated in the legend,
    # not in an in-panel text box that can obscure the axis.
    f.subplots_adjust(left=.078,right=.985,bottom=.09,top=.94)
    save(f,2,600)


def fig3():
    f=plt.figure(figsize=(7.0866,5.62)); outer=f.add_gridspec(2,1,height_ratios=[1.55,1],hspace=.56)
    top=outer[0].subgridspec(1,10,wspace=.68); bottom=outer[1].subgridspec(1,4,width_ratios=[.92,1.0,1.12,1.18],wspace=.72)
    A=load_fig(3,"A"); ax=f.add_subplot(top[0,:7]); panel(ax,"A","MAGMA gene-level association",x=-.07)
    for ch,g in A.groupby("chr_num"): ax.scatter(g.x,g.neglog10p,s=1.2,c=SCI if int(ch)%2 else "#7CAEB2",alpha=.75,rasterized=True)
    ax.axhline(-math.log10(A.bonferroni_threshold.iloc[0]),c=TX,ls="--",lw=.7); mh=A[A.mhc_flag.astype(bool)];
    if len(mh): ax.scatter(mh.x,mh.neglog10p,s=1.1,c="#BFC4C7",alpha=.45,rasterized=True)
    mids=A.groupby("chr_num").x.median(); ax.set_xticks(mids,mids.index.astype(int)); ax.set(xlabel="Chromosome",ylabel="-log10(P)"); tidy(ax,None)
    for tick in ax.get_xticklabels(): tick.set_fontsize(5.8)
    for tick in ax.get_xticklabels()[-3:]: tick.set_fontsize(4.2)
    labels=["DCC","EYS","FGFR3","SOX5","ANKS1B","PARK2","MAPK3","NFU1","BSN","GPR1","TXNL1","PDPR","GFPT1"]
    offsets={"NFU1":(-12,18),"GPR1":(8,26),"BSN":(-5,38),"FGFR3":(2,10),"EYS":(-5,15),"PARK2":(5,8),"SOX5":(-8,8),"ANKS1B":(8,14),"MAPK3":(-6,11),"PDPR":(-8,7),"TXNL1":(7,14),"GFPT1":(10,2),"DCC":(0,12)}
    for _,r in A[A.gene_symbol.isin(labels)].iterrows():
        dx,dy=offsets.get(r.gene_symbol,(0,7)); ax.annotate(r.gene_symbol,(r.x,r.neglog10p),xytext=(dx,dy),textcoords="offset points",ha="center",fontsize=6.0 if r.gene_symbol in GENE_COL else 6.1,c=GENE_COL.get(r.gene_symbol,DARK),fontweight="bold" if r.gene_symbol in GENE_COL else "normal",arrowprops=dict(arrowstyle="-",lw=.45,color=GENE_COL.get(r.gene_symbol,"#899399")))
    B=load_fig(3,"B"); B=B[~B.mhc_flag.astype(bool)].nsmallest(15,"p").sort_values("p",ascending=False)
    ax=f.add_subplot(top[0,7:]); panel(ax,"B","Top non-MHC genes",x=-.15)
    y=np.arange(len(B)); v=-np.log10(B.p); cols=[GENE_COL.get(g,"#8EA0A7") for g in B.gene_symbol]; ax.hlines(y,0,v,colors=LIGHT,lw=2); ax.scatter(v,y,c=cols,s=17,zorder=3); ax.set_yticks(y,B.gene_symbol,fontsize=5.7); ax.set_xlabel("-log10(P)"); tidy(ax,"x")
    P=pd.read_csv(ROOT/"results/phase3_final/03_magma/pathway_results.tsv",sep="\t",low_memory=False); P=P[(P.trait.eq("SCIATICA"))&(P.analysis.eq("primary"))&(P.type.eq("SET"))].copy()
    P["pnum"]=pd.to_numeric(P.p,errors="coerce"); P["qnum"]=pd.to_numeric(P.get("fdr"),errors="coerce"); P=P.dropna(subset=["pnum"]).sort_values("pnum").reset_index(drop=True); P["rank"]=np.arange(1,len(P)+1); P["nlp"]=-np.log10(P.pnum)
    ax=f.add_subplot(bottom[0,0]); panel(ax,"C","Pathway rank",x=-.18)
    ax.scatter(P["rank"],P.nlp,s=2,c="#B8C2C6",alpha=.5,rasterized=True); hit=P[P.qnum.lt(.05)]; ax.scatter(hit["rank"],hit.nlp,s=19,c=TX,zorder=4); ax.set(xlabel="Pathway rank",ylabel="-log10(P)"); tidy(ax)
    if len(hit):
        label=str(hit.iloc[0].get("full_name",hit.iloc[0].get("variable","FDR pathway"))).replace("_"," ")
        label="\n".join(textwrap.wrap(label,width=22))
        ax.annotate(label,(hit.iloc[0]["rank"],hit.iloc[0].nlp),xytext=(7,-1),textcoords="offset points",fontsize=6.0,c=TX,va="center",arrowprops=dict(arrowstyle="-",lw=.45,color=TX))
    D=load_fig(3,"D"); D["x"]=-np.log10(D.p_SCIATICA); D["y"]=-np.log10(D.p_LDH)
    D["class"]="Neither"; D.loc[D.fdr_SCIATICA.lt(.05)&~D.fdr_LDH.lt(.05),"class"]="Sciatica-first"; D.loc[D.fdr_SCIATICA.lt(.05)&D.fdr_LDH.lt(.05),"class"]="Shared"; D.loc[~D.fdr_SCIATICA.lt(.05)&D.fdr_LDH.lt(.05),"class"]="LDH-dominant"
    cc={"Neither":"#C8CED1","Sciatica-first":SCI,"Shared":"#6B6C78","LDH-dominant":LDH}; ax=f.add_subplot(bottom[0,1]); panel(ax,"D","Sciatica versus LDH",x=-.17)
    for c,g in D.groupby("class"): ax.scatter(g.x,g.y,s=2.4,c=cc[c],alpha=.55,label=c,rasterized=True)
    d_offsets={"TXNL1":(5,-7),"MAPK3":(5,11),"FGFR3":(-8,7),"PDPR":(8,-8),"GFPT1":(-8,-7)}
    for _,r in D[D.gene_symbol.isin(GENE_COL)].iterrows():
        dx,dy=d_offsets[r.gene_symbol]; ax.annotate(r.gene_symbol,(r.x,r.y),xytext=(dx,dy),textcoords="offset points",fontsize=6.0,c=GENE_COL[r.gene_symbol],fontweight="bold",arrowprops=dict(arrowstyle="-",lw=.4,color=GENE_COL[r.gene_symbol]))
    ax.set(xlabel="Sciatica -log10(P)",ylabel="LDH -log10(P)"); tidy(ax,None)
    ax=f.add_subplot(bottom[0,2]); panel(ax,"E","Gene-class composition",x=-.17)
    cnt=D["class"].value_counts().reindex(["Sciatica-first","Shared","LDH-dominant","Neither"]).fillna(0); ax.barh(np.arange(4),cnt,color=[SCI,"#6B6C78",LDH,"#CDD3D6"]); ax.set_yticks(np.arange(4),cnt.index,fontsize=5.8); ax.invert_yaxis(); ax.set_xlabel("Genes"); tidy(ax,"x")
    E=load_fig(3,"E").head(12).sort_values("p_primary",ascending=False); ax=f.add_subplot(bottom[0,3]); panel(ax,"F","MHC-excluded\nsensitivity",x=-.17,title_pad=8)
    y=np.arange(len(E)); x1=-np.log10(E.p_primary); x2=-np.log10(E.p_no_mhc); ax.hlines(y,x1,x2,colors=LIGHT,lw=1); ax.scatter(x1,y,c=SCI,s=12,label="Primary"); ax.scatter(x2,y,c=TX,s=12,label="No MHC"); ax.set_yticks(y,E.gene_symbol,fontsize=6.4); ax.tick_params(axis="y",pad=2); ax.set_xlabel("-log10(P)"); ax.legend(frameon=False,loc="lower right",fontsize=6.2); tidy(ax,"x")
    for name,d in [("A_MAGMA",A),("B_top_genes",B),("C_pathways",P),("D_E_gene_classes",D),("F_MHC_sensitivity",E)]: write_source(3,name,d)
    f.subplots_adjust(left=.085,right=.985,bottom=.095,top=.94)
    save(f,3,600)


def fig4():
    f=plt.figure(figsize=(7.0866,5.78)); outer=f.add_gridspec(3,12,height_ratios=[1.35,.9,1.25],hspace=.72,wspace=.92)
    A=load_fig(4,"A"); ax=f.add_subplot(outer[0,:8]); panel(ax,"A","Five-tissue FUSION landscape",x=-.07)
    for t,g in A.groupby("tissue"): ax.scatter(g.x,g.neglog10p,s=.8,c=TISSUE_COL.get(t,GREY),alpha=.65,label=t.replace("_"," "),rasterized=True)
    ax.axhline(-math.log10(.05/len(A)),c=TX,ls="--",lw=.65); mids=A.groupby("CHR").x.median(); ax.set_xticks(mids,mids.index.astype(int)); ax.set(xlabel="Chromosome",ylabel="-log10(P)"); tidy(ax,None)
    for tick in ax.get_xticklabels(): tick.set_fontsize(5.8)
    for tick in ax.get_xticklabels()[-3:]: tick.set_fontsize(4.2)
    handles=[Line2D([0],[0],marker="o",ls="",ms=3.5,color=TISSUE_COL[t],label=TISSUE_LABEL[t]) for t in TISSUE_COL]
    ax.legend(handles=handles,ncol=3,frameon=False,loc="upper left",bbox_to_anchor=(0,.99),handletextpad=.25,columnspacing=.7,fontsize=6.0,borderaxespad=0)
    cand=["TXNL1","MAPK3","FGFR3","PDPR","GFPT1"]; H=A[A.gene_name.isin(cand)].pivot_table(index="gene_name",columns="tissue",values="zscore",aggfunc="first").reindex(cand)
    ax=f.add_subplot(outer[0,8:]); panel(ax,"B","Candidate × tissue Z",x=-.15)
    arr=H.reindex(columns=TISSUE_COL).to_numpy(); im=ax.imshow(np.ma.masked_invalid(arr),cmap="RdBu_r",norm=TwoSlopeNorm(0,vmin=np.nanmin(arr),vmax=np.nanmax(arr)),aspect="auto"); ax.set_yticks(range(5),cand); ax.set_xticks(range(5),[TISSUE_LABEL[x] for x in TISSUE_COL],rotation=45,ha="right",fontsize=6.1)
    for i in range(5):
        for j in range(5):
            if np.isnan(arr[i,j]): ax.text(j,i,"×",ha="center",va="center",c=GREY,fontsize=7)
            else: ax.text(j,i,f"{arr[i,j]:.1f}",ha="center",va="center",fontsize=4.8,c="white" if abs(arr[i,j])>4 else DARK)
    B=load_fig(4,"B"); top=B[~B.one_snp_model.astype(bool)].nsmallest(12,"pvalue").sort_values("zscore")
    ax=f.add_subplot(outer[1,:6]); panel(ax,"C","Top multi-SNP models",x=-.08)
    short=lambda t: TISSUE_LABEL.get(t,t.replace("_"," "))
    y=np.arange(len(top)); cols=[GENE_COL.get(g,TISSUE_COL.get(t,"#8AA0A6")) for g,t in zip(top.gene_name,top.tissue)]; ax.hlines(y,0,top.zscore,colors=LIGHT,lw=2); ax.scatter(top.zscore,y,c=cols,s=18); ax.axvline(0,c=GREY,lw=.55); ax.set_yticks(y,top.gene_name,fontsize=6.2); ax.tick_params(axis="y",pad=2); ax.set_xlabel("FUSION Z"); tidy(ax,"x")
    C=load_fig(4,"C").sort_values("JOINT.Z"); ax=f.add_subplot(outer[1,6:]); panel(ax,"D","Eleven joint-retained models",x=-.08)
    y=np.arange(len(C)); ax.hlines(y,0,C["JOINT.Z"],colors=LIGHT,lw=2); ax.scatter(C["JOINT.Z"],y,c=[GENE_COL.get(g,SCI) for g in C.gene_name],s=19); ax.axvline(0,c=GREY,lw=.55); ax.set_yticks(y,C.gene_name,fontsize=5.4); ax.set_xlabel("Joint Z"); tidy(ax,"x")
    bottom=outer[2,:].subgridspec(1,12,wspace=1.35)
    L=load_fig(4,"D").sort_values("position_gwas"); ax=f.add_subplot(bottom[0,:6]); panel(ax,"E","TXNL1 full-cis locus",x=-.08,title_pad=10)
    ax2=ax.twinx(); ax.scatter(L.position_gwas/1e6,L.gwas_nlp,s=2.5,c=SCI,alpha=.65,rasterized=True,label="GWAS"); ax2.scatter(L.position_qtl/1e6,L.qtl_nlp,s=2.5,c=TX,alpha=.45,rasterized=True,label="BLUEPRINT eQTL"); ax.set(xlabel="chr18 position (Mb)",ylabel="GWAS -log10(P)"); ax2.tick_params(colors=TX,labelsize=4.6,pad=-7,direction="in"); ax2.set_ylabel(""); tidy(ax,None); ax2.spines["right"].set_visible(True)
    ax.tick_params(axis="both",labelsize=5.4,pad=2); ax.yaxis.label.set_size(6.1); ax.text(.98,.94,"ABF PP4 = 0.888\nSuSiE PP4 = 0.973",transform=ax.transAxes,ha="right",va="top",fontsize=5.7); ax.text(.98,.04,"LD/credible set unavailable",transform=ax.transAxes,ha="right",va="bottom",fontsize=5.4,color=GREY)
    E=load_fig(4,"E").set_index("gene").reindex(cand).reset_index(); ax=f.add_subplot(bottom[0,7:9]); panel(ax,"F","ABF versus\nSuSiE",x=-.18,title_pad=7)
    yy=np.arange(len(E));
    for i,r in E.iterrows():
        if pd.notna(r.SuSiE_max_PP4): ax.plot([r.ABF_max_PP4,r.SuSiE_max_PP4],[i,i],c=LIGHT,lw=2)
    ax.scatter(E.ABF_max_PP4,yy,c=SCI,s=18,label="ABF"); ok=E.SuSiE_max_PP4.notna(); ax.scatter(E.loc[ok,"SuSiE_max_PP4"],yy[ok],c=TX,s=18,label="SuSiE"); ax.axvline(.8,c=GREY,ls="--",lw=.6); ax.set_yticks(yy,E.gene,fontsize=5.5); ax.tick_params(axis="y",pad=1.2,labelsize=5.5); ax.tick_params(axis="x",labelsize=5.4); ax.set(xlim=(.75,1.01),xlabel="PP4"); ax.invert_yaxis(); ax.legend(frameon=False,loc="lower right",fontsize=5.5); tidy(ax,"x")
    cov=load_fig(4,"G"); ax=f.add_subplot(bottom[0,9:]); panel(ax,"G","Context coverage",x=-.15,y=1.10,title_pad=10)
    ax.bar(np.arange(5)-.17,cov.tested,.34,color=LIGHT,label="Usable"); ax.bar(np.arange(5)+.17,cov.strong,.34,color=SCI,label="Strong"); ax.set_xticks(range(5),cov.gene,rotation=45,ha="right",fontsize=5.6); ax.set_ylabel("Contexts",labelpad=1); ax.set_ylim(0,max(cov.tested.max()+3.2,10)); ax.legend(frameon=False,fontsize=5.4,ncol=2,loc="upper center",bbox_to_anchor=(.52,.94),columnspacing=.45,handletextpad=.22,borderaxespad=0); tidy(ax)
    for name,d in [("A_FUSION",A),("B_candidate_heatmap",H.reset_index()),("C_top_multisnp",top),("D_joint",C),("E_TXNL1_locus",L),("F_coloc",E),("G_context_coverage",cov)]: write_source(4,name,d)
    f.subplots_adjust(left=.105,right=.985,bottom=.105,top=.94)
    save(f,4,600)


def fig5():
    f=plt.figure(figsize=(7.0866,6.05)); gs=f.add_gridspec(3,12,height_ratios=[1.45,1.0,1.0],hspace=.88,wspace=1.02)
    maps=[load_fig(5,x) for x in "BCDE"]; allz=pd.concat(maps).z; vmax=float(np.nanquantile(abs(allz),.99)); map_axes=[]
    for j,(d,letter) in enumerate(zip(maps,"ABCD")):
        ax=f.add_subplot(gs[0,j*3:(j+1)*3]); map_axes.append(ax); panel(ax,letter,None,x=-.13,y=1.18)
        ax.text(0,1.14,str(d.donor.iloc[0]),transform=ax.transAxes,ha="left",va="bottom",fontsize=7.6,fontweight="semibold",color=DARK,clip_on=False)
        ax.text(0,1.03,str(d['section'].iloc[0]),transform=ax.transAxes,ha="left",va="bottom",fontsize=6.2,color=GREY,clip_on=False)
        sc=ax.scatter(d.spatial_x,d.spatial_y,c=d.z,s=2.6,cmap="RdBu_r",vmin=-vmax,vmax=vmax,rasterized=True); ax.set_aspect("equal"); ax.axis("off")
    cb=f.colorbar(sc,ax=map_axes,orientation="horizontal",fraction=.035,pad=.055,aspect=55); cb.set_label("gsMap z",fontsize=6.5,labelpad=2); cb.ax.tick_params(labelsize=6.2,pad=1)
    D=load_fig(5,"F"); ax=f.add_subplot(gs[1,:3]); panel(ax,"E","Donor-level evidence",x=-.18)
    ax.hlines(np.arange(4),0,-np.log10(D.p_acat_sections),colors=LIGHT,lw=2); ax.scatter(-np.log10(D.p_acat_sections),np.arange(4),c=SCI,s=24); ax.axvline(-math.log10(.05),c=GREY,ls="--",lw=.6); ax.set_yticks(range(4),D.donor); ax.set_xlabel("-log10(P)"); tidy(ax,"x")
    methods=load_fig(5,"G"); methods=methods[~methods.method.astype(str).str.contains("Stouffer",case=False,na=False)].copy(); ax=f.add_subplot(gs[1,3:6]); panel(ax,"F","Cross-donor aggregation",x=-.18)
    valid=methods.p.notna(); ax.bar(np.arange(len(methods))[valid],-np.log10(methods.loc[valid,"p"]),color=[SCI,"#72A59E"][:int(valid.sum())])
    ax.set_xticks(range(len(methods)),methods.method.str.replace(" (primary)","",regex=False),rotation=25,ha="right",fontsize=6.0); ax.set_ylabel("-log10(P)"); tidy(ax)
    H=load_fig(5,"H"); ax=f.add_subplot(gs[1,6:9]); panel(ax,"G","Leave-one-donor-out",x=-.18)
    ax.scatter(np.arange(4),-np.log10(H.p_acat),c=SCI,s=25); ax.axhline(-math.log10(.05),c=GREY,ls="--",lw=.6); ax.set_xticks(range(4),H.left_out_donor.str.replace("Donor","D")); ax.set_ylabel("-log10(P)"); tidy(ax)
    I=load_fig(5,"I"); di=I[I.level.astype(str).str.contains("donor",case=False,na=False)]; ax=f.add_subplot(gs[1,9:]); panel(ax,"H","Sciatica–LDH concordance",x=-.18)
    ax.scatter(np.arange(len(di)),di.spearman_z,c=LDH,s=21); med=di.spearman_z.median(); ax.axhline(med,c=DARK,lw=.7); ax.set_xticks(np.arange(len(di)),di.donor.str.replace("Donor","D"),rotation=45); ax.set_ylabel("Spearman rho"); ax.text(.97,.93,f"median {med:.2f}",transform=ax.transAxes,ha="right",va="top",fontsize=6.0,color=GREY); tidy(ax)
    lower=gs[2,:].subgridspec(1,12,wspace=1.20)
    J=load_fig(5,"J"); ax=f.add_subplot(lower[0,:4]); panel(ax,"I","Sensitivity and calibration",x=-.18)
    vals=[]
    if "max_abs_z_delta" in J: vals.append(("MHC\nexclusion",pd.to_numeric(J.max_abs_z_delta,errors="coerce").max()))
    hp=pd.to_numeric(J.get("global_p_cauchy"),errors="coerce"); vals.append(("Height trait\n−log10(P)",-math.log10(hp.dropna().iloc[0]) if hp.notna().any() and hp.dropna().iloc[0]>0 else np.nan)); vals.append(("Common variants\nnot available","NA"))
    for i,(lab,val) in enumerate(vals):
        ax.add_patch(FancyBboxPatch((i-.40,.12),.80,.64,boxstyle="round,pad=.025",fc=PALE,ec=LIGHT)); ax.text(i,.55,lab,ha="center",fontsize=5.7); ax.text(i,.29,f"{val:.3g}" if isinstance(val,float) and np.isfinite(val) else "NA",ha="center",fontsize=5.9,fontweight="bold",c=SCI if i<2 else GREY)
    ax.set(xlim=(-.6,2.6),ylim=(0,1)); ax.axis("off")
    k=load_fig(5,"K"); ax=f.add_subplot(lower[0,5:12]); panel(ax,"J","TXNL1 spatial coherence across 20 sections",x=-.12)
    donors=sorted(k.donor.unique()); shapes=["o","s","^","D"]
    for do,sh in zip(donors,shapes):
        q=k[k.donor.eq(do)]; ax.scatter(q.fraction_expressing,q.spearman_rho,s=20,marker=sh,label=do,c=SCI,alpha=.8)
    ax.axhline(0,c=GREY,lw=.55); ax.set(xlabel="TXNL1 fraction detected",ylabel="Section-level Spearman rho"); ax.xaxis.label.set_size(6.2); ax.yaxis.label.set_size(6.2); ax.tick_params(labelsize=5.8); ax.legend(frameon=False,ncol=2,loc="upper center",bbox_to_anchor=(.52,-.28),fontsize=5.7,borderaxespad=0,columnspacing=.7,handletextpad=.3); tidy(ax)
    # Median/recurrence statistics remain in the legend and source data rather than a box over the scatter.
    J_source=J.copy()
    for col in J_source.select_dtypes(include="object").columns:
        J_source[col]=J_source[col].astype(str).str.replace("height LD-preserving control","height LD-preserving calibration trait",case=False,regex=False).str.replace("negative control","calibration trait",case=False,regex=False)
    for name,d in [("B_E_maps",pd.concat(maps)),("F_donor",D),("G_meta",methods),("H_LODO",H),("I_concordance",di),("J_controls",J_source),("K_TXNL1",k)]: write_source(5,name,d)
    f.subplots_adjust(left=.09,right=.975,bottom=.15,top=.93)
    save(f,5,300)


def fig6():
    f=plt.figure(figsize=(7.0866,6.32)); gs=f.add_gridspec(3,14,height_ratios=[1.48,1.0,.92],hspace=.94,wspace=1.12)
    U=load_fig(6,"A"); ax=f.add_subplot(gs[:2,:5]); panel(ax,"A",None,x=-.08); ax.text(0,1.035,"Adult spinal-cord snRNA atlas",transform=ax.transAxes,ha="left",va="bottom",fontsize=7.8,fontweight="semibold",color=DARK,clip_on=False)
    direct={"Inhibitory neuron","Excitatory neuron","Astrocyte","Microglia","Oligodendrocyte"}; minor=[]
    for c,g in U.groupby("cell_type"):
        ax.scatter(g.UMAP1,g.UMAP2,s=.34,c=CELL_COL.get(c,"#B8C0C4"),alpha=.65,rasterized=True)
        if c in direct:
            ax.text(g.UMAP1.median(),g.UMAP2.median(),str(c).replace(" neuron","\nneuron"),fontsize=6.0,ha="center",va="center",c=DARK,bbox=dict(fc="white",ec="none",alpha=.68,pad=.8))
        else:
            minor.append(c)
    ax.set(xticks=[],yticks=[],xlabel="UMAP1",ylabel="UMAP2"); ax.text(.02,.02,"68,175 nuclei · 9 donors",transform=ax.transAxes,fontsize=6.3,fontweight="bold")
    if minor:
        handles=[Line2D([0],[0],marker="o",ls="",ms=3.5,color=CELL_COL.get(c,"#B8C0C4"),label=CELL_SHORT.get(c,c)) for c in minor]
        ax.legend(handles=handles,ncol=2,frameon=True,facecolor="white",edgecolor=LIGHT,loc="upper left",fontsize=6.0,handletextpad=.25,columnspacing=.5,borderpad=.3)
    markers=load_fig(6,"markers"); major=list(U.cell_type.value_counts().head(8).index); mm=markers[markers.cell_type.isin(major)&markers.specificity_rank.le(4)].copy(); mm["score"]=5-mm.specificity_rank
    piv=mm.pivot_table(index="cell_type",columns="gene",values="score",fill_value=0).reindex([x for x in major if x in mm.cell_type.unique()]); piv=piv.loc[:,piv.max().sort_values(ascending=False).index[:24]]
    ax=f.add_subplot(gs[0,6:10]); panel(ax,"B",None,x=-.12); ax.text(0,1.035,"Marker audit",transform=ax.transAxes,ha="left",va="bottom",fontsize=7.8,fontweight="semibold",color=DARK,clip_on=False)
    marker_rows=[CELL_SHORT.get(str(x),str(x)) for x in piv.index]
    ax.imshow(piv,cmap="Blues",aspect="auto",vmin=0,vmax=4); ax.set_yticks(range(len(piv)),marker_rows,fontsize=5.3); ax.tick_params(axis="y",pad=1); ax.set_xticks(range(len(piv.columns)),piv.columns,rotation=70,ha="right",fontsize=5.7); ax.tick_params(axis="x",pad=2)
    B=load_fig(6,"B"); ax=f.add_subplot(gs[0,11:]); panel(ax,"C",None,x=-.18); ax.text(0,1.035,"Donor composition",transform=ax.transAxes,ha="left",va="bottom",fontsize=7.8,fontweight="semibold",color=DARK,clip_on=False)
    top_cells=list(B.groupby("cell_type").fraction.sum().sort_values(ascending=False).head(7).index)
    bc=B.copy(); bc["display_cell_type"]=bc.cell_type.where(bc.cell_type.isin(top_cells),"Other")
    pp=bc.pivot_table(index="donor",columns="display_cell_type",values="fraction",aggfunc="sum",fill_value=0); bottom=np.zeros(len(pp))
    for c in pp.columns:
        ax.bar(range(len(pp)),pp[c],bottom=bottom,color=CELL_COL.get(c,"#B8C0C4"),width=.8,label=CELL_SHORT.get(c,c)); bottom+=pp[c].values
    ax.set_xticks(range(len(pp)),[str(x).replace("Donor","D") for x in pp.index],rotation=45,ha="right",fontsize=6.0); ax.set_ylabel("Fraction"); ax.legend(frameon=False,bbox_to_anchor=(1.0,1),loc="upper left",fontsize=5.8,ncol=1,columnspacing=.35,handletextpad=.25,borderaxespad=0); tidy(ax,None)
    allassoc=load_fig(6,"assoc"); order=allassoc.groupby("cell_type").spearman_rho.median().sort_values().index
    ax=f.add_subplot(gs[1,6:]); panel(ax,"D",None,x=-.10); ax.text(0,1.035,"gsMap–cell-type association across sections",transform=ax.transAxes,ha="left",va="bottom",fontsize=7.8,fontweight="semibold",color=DARK,clip_on=False)
    data=[allassoc.loc[allassoc.cell_type.eq(c),"spearman_rho"] for c in order]; bp=ax.boxplot(data,vert=False,positions=np.arange(len(order)),widths=.5,patch_artist=True,showfliers=False,medianprops=dict(color=DARK,lw=1));
    for patch,c in zip(bp["boxes"],order): patch.set_facecolor(CELL_COL.get(c,LIGHT)); patch.set_alpha(.65)
    xmin=min(-.45,float(allassoc.spearman_rho.min())-.03); xmax=max(.35,float(allassoc.spearman_rho.max())+.10); ax.set_xlim(xmin,xmax)
    count_x=xmax-.018
    for i,c in enumerate(order):
        q=allassoc[allassoc.cell_type.eq(c)]; ax.scatter(q.spearman_rho,np.full(len(q),i),s=5,c=CELL_COL.get(c,GREY),alpha=.5)
        ax.text(count_x,i,f"{(q.spearman_rho>0).sum()}/20",fontsize=6.0,va="center",ha="right",bbox=dict(fc="white",ec="none",alpha=.7,pad=.3))
    assoc_rows=[CELL_SHORT.get(str(x),str(x)) for x in order]
    ax.axvline(0,c=GREY,lw=.6); ax.set_yticks(np.arange(len(order)),assoc_rows,fontsize=5.3); ax.tick_params(axis="y",pad=1); ax.set_xlabel("Section-level Spearman rho"); tidy(ax,"x")
    cand=["TXNL1","MAPK3","FGFR3","PDPR","GFPT1"]; C=load_fig(6,"C"); E=load_fig(6,"E")
    bottom=gs[2,:].subgridspec(1,100,wspace=2.8)
    for loc,dat,title,letter in [(bottom[0,:18],C,"Spinal candidate\nexpression","E"),(bottom[0,21:35],E,"DRG candidate\nexpression","F")]:
        ax=f.add_subplot(loc); panel(ax,letter,None,x=-.16 if letter=="F" else -.13); ax.text(0,1.035,title,transform=ax.transAxes,ha="left",va="bottom",fontsize=7.3,fontweight="semibold",color=DARK,clip_on=False); q=dat[dat.gene.isin(cand)].copy(); genes=[g for g in cand if g in q.gene.unique()]; cells=list(q.groupby("cell_type").fraction_expressing.mean().sort_values(ascending=False).head(7).index); q=q[q.cell_type.isin(cells)];
        for _,r in q.iterrows(): ax.scatter(cells.index(r.cell_type),genes.index(r.gene),s=max(4,r.fraction_expressing*145),c=r.mean_log_expression,cmap="viridis",vmin=q.mean_log_expression.min(),vmax=q.mean_log_expression.max())
        ax.set_xticks(range(len(cells)),[x.replace(" ","\n") for x in cells],rotation=60 if letter=="E" else 45,ha="right",fontsize=6.0); ax.set_yticks(range(len(genes)),[]); ax.set_xlim(-1.75,len(cells)-.5); ax.invert_yaxis(); tidy(ax,None)
        for yi,gene in enumerate(genes): ax.text(-1.65,yi,gene,ha="left",va="center",fontsize=6.4,color=DARK)
        if letter=="F":
            ax.set_ylim(len(genes)-.5,-1.0)
            ax.text(-1.65,-.72,"6 prep. / 5 unique donors",ha="left",va="center",fontsize=6.0,c=GREY)
    reg=load_fig(6,"registry").set_index("gene").reindex(cand).reset_index()
    recurrence=pd.read_csv(DRG_RECURRENCE,sep="\t")
    sensory=recurrence[recurrence.cell_type.isin(["Nociceptor","Mechanoreceptor"])].copy()
    drg_rec=(sensory.groupby("gene",as_index=False)
             .agg(DRG_positive_unique_donors=("expression_positive_unique_donors","max"),
                  DRG_assessable_unique_donors=("assessable_unique_donors","max"),
                  DRG_positive_preparations=("expression_positive_preparations","max"),
                  DRG_assessable_preparations=("assessable_preparations","max")))
    reg=reg.merge(drg_rec,on="gene",how="left")
    ax=f.add_subplot(bottom[0,38:52]); panel(ax,"G",None,x=-.22); ax.text(0,1.035,"Donor recurrence",transform=ax.transAxes,ha="left",va="bottom",fontsize=6.7,fontweight="semibold",color=DARK,clip_on=False)
    y=np.arange(3); ax.barh(y-.18,reg.loc[:2,"snRNA_spinal_donors"],height=.32,color=SCI,label="Spinal"); ax.barh(y+.18,reg.loc[:2,"DRG_positive_unique_donors"],height=.32,color=TX,label="DRG sensory"); ax.set_yticks(y,[]); ax.set_xlim(-3.0,10); ax.set_ylim(2.5,-1.0); ax.legend(frameon=False,fontsize=6.0,ncol=2,loc="upper center",bbox_to_anchor=(.62,.99),columnspacing=.6,handletextpad=.3); tidy(ax,"x")
    for yi,gene in enumerate(reg.loc[:2,"gene"]): ax.text(-2.85,yi,gene,ha="left",va="center",fontsize=6.3,color=DARK)
    for yi,r in reg.loc[:2].iterrows():
        ax.text(r.snRNA_spinal_donors+.15,yi-.18,f"{int(r.snRNA_spinal_donors)}/9",fontsize=5.4,va="center",color=SCI)
        ax.text(r.DRG_positive_unique_donors+.15,yi+.18,f"{int(r.DRG_positive_unique_donors)}/{int(r.DRG_assessable_unique_donors)}",fontsize=5.4,va="center",color=TX)
    ax=f.add_subplot(bottom[0,56:]); panel(ax,"H",None,x=-.12); ax.text(0,1.035,"Categorical evidence matrix",transform=ax.transAxes,ha="left",va="bottom",fontsize=6.9,fontweight="semibold",color=DARK,clip_on=False)
    show=reg.loc[:,["gene","SCIATICA_MAGMA_FDR","robust_TWAS_z","coloc_PP4","susie_coloc_max_PP4","gsMap_positive_sections","snRNA_spinal_donors","DRG_positive_unique_donors","DRG_assessable_unique_donors"]].copy(); mat=np.zeros((len(show),7)); mat[:,0]=show.SCIATICA_MAGMA_FDR.lt(.05); mat[:,1]=show.robust_TWAS_z.notna(); mat[:,2]=show.coloc_PP4.ge(.8); mat[:,3]=show.susie_coloc_max_PP4.ge(.8); mat[:,4]=show.gsMap_positive_sections/20; mat[:,5]=show.snRNA_spinal_donors/9; mat[:,6]=show.DRG_positive_unique_donors/show.DRG_assessable_unique_donors
    ax.imshow(mat,cmap=LinearSegmentedColormap.from_list("cat",["#E5E8EA","#7CB8AD",SCI]),vmin=0,vmax=1,aspect="auto"); ax.set_xticks(range(7),["MAGMA","TWAS","ABF","SuSiE","Spatial","Spinal","DRG"],rotation=55,ha="right",fontsize=6.0); ax.set_yticks(range(5),[]); ax.set_xlim(-2.8,6.5)
    for yi,gene in enumerate(show.gene): ax.text(-2.7,yi,gene,ha="left",va="center",fontsize=7.0,color=DARK)
    for i,r in show.iterrows():
        txt=[f"{r.SCIATICA_MAGMA_FDR:.1e}".replace("e-","e−"),f"{r.robust_TWAS_z:.1f}".replace("-","−"),f"{r.coloc_PP4:.2f}","NA" if pd.isna(r.susie_coloc_max_PP4) else f"{r.susie_coloc_max_PP4:.2f}",f"{int(r.gsMap_positive_sections)}/20",f"{int(r.snRNA_spinal_donors)}/9",("NA" if pd.isna(r.DRG_assessable_unique_donors) or r.DRG_assessable_unique_donors==0 else f"{int(r.DRG_positive_unique_donors)}/{int(r.DRG_assessable_unique_donors)}")]
        for j,t in enumerate(txt): ax.text(j,i,t,ha="center",va="center",fontsize=5.35,c="white" if mat[i,j]>.75 else DARK)
    ax.add_patch(FancyBboxPatch((0,-.86),1,.30,boxstyle="round,pad=.015",fc=PALE,ec=LIGHT,transform=ax.transAxes,clip_on=False))
    ax.text(.02,-.71,"Final hierarchy",transform=ax.transAxes,fontsize=6.7,fontweight="bold",c=NAVY,clip_on=False,va="center")
    ax.text(.38,-.71,"1 TXNL1 — primary\n2 MAPK3 — backup   3 FGFR3 — third",transform=ax.transAxes,fontsize=6.2,c=DARK,clip_on=False,va="center")
    for name,d in [("A_UMAP",U),("B_markers",mm),("C_composition",B),("D_celltype_association",allassoc),("E_spinal_expression",C),("F_DRG_expression",E),("G_H_registry",reg)]: write_source(6,name,d)
    shrink_figure_elements(f,text_factor=.84,marker_factor=.78,min_text=5.25)
    f.subplots_adjust(left=.105,right=.845,bottom=.255,top=.91)
    save(f,6,600)




if __name__ == "__main__":
    setup()
    for fn in (fig1,fig2,fig3,fig4,fig5,fig6):
        print("building", fn.__name__, flush=True)
        fn()
    print(OUT)
