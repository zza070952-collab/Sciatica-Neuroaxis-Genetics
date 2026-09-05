from __future__ import annotations

import math
import os
import shutil
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
import fitz

ROOT = Path(os.environ.get("SCIATICA_PROJECT_ROOT", Path.cwd()))
PHASE = Path(os.environ.get("PHASE8_OUT", ROOT / "results/phase8_JTM_core_final"))
SRC = Path(os.environ.get("PHASE8_SOURCE_DATA", ROOT / "data/computational_source_data"))
OUT = PHASE / "07_supplementary_figures"
SD = OUT / "Supplementary_Figure_Source_Data"
OUT.mkdir(parents=True, exist_ok=True)
SD.mkdir(parents=True, exist_ok=True)

SCI = "#1F6F78"
LDH = "#B9852B"
TX = "#B65346"
BLUE = "#356A9A"
GREEN = "#668B63"
GREY = "#697680"
LIGHT = "#E5EBEE"
DARK = "#263740"


def style():
    mpl.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 6.8,
        "axes.titlesize": 7.8, "axes.labelsize": 6.8,
        "xtick.labelsize": 6.0, "ytick.labelsize": 6.0,
        "legend.fontsize": 5.8, "axes.spines.top": False,
        "axes.spines.right": False, "pdf.fonttype": 42,
        "ps.fonttype": 42, "svg.fonttype": "none",
    })


def panel(ax, letter, title):
    ax.text(-0.12, 1.04, letter, transform=ax.transAxes, fontsize=9.5,
            fontweight="bold", color=DARK, va="bottom")
    ax.set_title(title, loc="left", fontweight="bold", pad=7)


def tidy(ax, grid="y"):
    if grid:
        ax.grid(axis=grid, color="#DDE4E7", lw=0.45, zorder=0)
    ax.tick_params(length=2)


def load_fig(n, sheet):
    return pd.read_excel(SRC / f"Source_Data_Figure{n}.xlsx", sheet_name=sheet)


def save(fig, num):
    fig.subplots_adjust(left=.075, right=.98, bottom=.08, top=.95, hspace=.62, wspace=.55)
    base = OUT / f"Supplementary_Figure_S{num}"
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight", pad_inches=.035)
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight", pad_inches=.035)
    fig.savefig(base.with_suffix(".tiff"), dpi=600, bbox_inches="tight", pad_inches=.035,
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    return base.with_suffix(".pdf")


def s6():
    F = load_fig(5, "F_donor")
    G = load_fig(5, "G_meta")
    H = load_fig(5, "H_LODO")
    I = load_fig(5, "I_concordance")
    J = load_fig(5, "J_controls")
    G = G.loc[G["method"].astype(str).str.startswith("ACAT") | G["method"].astype(str).eq("Fisher")].copy()
    fig, aa = plt.subplots(2, 4, figsize=(7.0866, 5.6))
    ax = aa[0, 0]; panel(ax, "A", "Donor-level aggregation")
    ax.barh(F.donor, -np.log10(F.p_acat_sections), color=SCI); ax.axvline(-math.log10(.05), c=GREY, ls="--", lw=.7)
    ax.set_xlabel("−log10(P)"); tidy(ax, "x")
    ax = aa[0, 1]; panel(ax, "B", "Across-donor evidence")
    ax.bar(np.arange(len(G)), -np.log10(G.p), color=[SCI, "#74A99F"][:len(G)])
    ax.set_xticks(np.arange(len(G)), G.method); ax.set_ylabel("−log10(P)"); tidy(ax)
    ax = aa[0, 2]; panel(ax, "C", "Leave-one-donor-out")
    ax.scatter(np.arange(len(H)), -np.log10(H.p_acat), c=SCI, s=22)
    ax.axhline(-math.log10(.05), c=GREY, ls="--", lw=.7)
    ax.set_xticks(np.arange(len(H)), H.left_out_donor.astype(str).str.replace("Donor", "D")); ax.set_ylabel("−log10(P)"); tidy(ax)
    ax = aa[0, 3]; panel(ax, "D", "Sciatica–LDH concordance")
    z = I[I.level.astype(str).str.contains("donor", case=False, na=False)].copy()
    ax.scatter(z.spearman_z, z.top5_jaccard, c=LDH, s=24)
    for _, r in z.iterrows():
        ax.annotate(str(r.donor).replace("Donor", "D"), (r.spearman_z, r.top5_jaccard), xytext=(3, 3), textcoords="offset points", fontsize=5.7)
    ax.set(xlabel="Spearman ρ", ylabel="Top-5% Jaccard"); tidy(ax)

    ax = aa[1, 0]; panel(ax, "E", "MHC exclusion")
    mz = pd.to_numeric(J.get("max_abs_z_delta"), errors="coerce").max()
    ax.axis("off"); ax.text(.04, .78, f"Maximum |Δz|\n{mz:.3g}", fontsize=9, color=SCI, fontweight="bold", va="top")
    ax.text(.04, .35, "Sensitivity analysis;\nthe main spatial conclusion\nwas unchanged.", fontsize=6.0, va="top", linespacing=1.35)
    ax = aa[1, 1]; panel(ax, "F", "Height calibration trait")
    hp = pd.to_numeric(J.get("global_p_cauchy"), errors="coerce").dropna()
    hv = -math.log10(hp.iloc[0]) if len(hp) and hp.iloc[0] > 0 else np.nan
    ax.bar([0], [hv], color=GREY, width=.5); ax.set_xticks([0], ["Height"]); ax.set_ylabel("−log10(P)"); tidy(ax)
    ax.text(.5, -.28, "LD-preserving reference trait;\nnot a null/negative control", transform=ax.transAxes, ha="center", fontsize=5.5)
    ax = aa[1, 2]; panel(ax, "G", "Replicate hierarchy")
    ax.axis("off"); ax.text(.05, .82, "20 sections\n4 adult donors\n85,100 spots", va="top", fontsize=10, color=SCI, linespacing=1.5)
    ax.text(.05, .25, "Donor was the biological replicate;\nspots were not treated as\nindependent replicates.", fontsize=5.9, va="top", linespacing=1.3)
    ax = aa[1, 3]; panel(ax, "H", "Interpretation boundary")
    ax.axis("off"); ax.text(.05, .83, "Normal adult reference atlas\nNo case–control expression test\nNo dorsal/ventral horn labels\nLDH was a structural reference", va="top", fontsize=7.1, linespacing=1.55)

    F.to_csv(SD / "S6_A_donor.tsv", sep="\t", index=False)
    G.to_csv(SD / "S6_B_meta.tsv", sep="\t", index=False)
    H.to_csv(SD / "S6_C_LODO.tsv", sep="\t", index=False)
    z.to_csv(SD / "S6_D_concordance.tsv", sep="\t", index=False)
    J.to_csv(SD / "S6_EF_sensitivity.tsv", sep="\t", index=False)
    return save(fig, 6)


def s8():
    expr = load_fig(6, "F_DRG_expression")
    rec = pd.read_csv(PHASE / "03_drg_recurrence_fix/DRG_expression_positive_recurrence.tsv", sep="\t")
    sens = pd.read_csv(PHASE / "03_drg_recurrence_fix/DRG_recurrence_threshold_sensitivity.tsv", sep="\t")
    genes = ["TXNL1", "MAPK3", "FGFR3", "PDPR", "GFPT1"]
    classes = ["nociceptor", "mechanoreceptor", "Schwann cell"]
    fig = plt.figure(figsize=(7.0866, 5.8))
    gs = fig.add_gridspec(2, 3, height_ratios=[.72, 1.45], width_ratios=[.92, 1.02, 1.36], wspace=.42)
    texts = [
        ("A", "Dataset coverage", "1,837 QC-passed nuclei\n6 preparations\n5 unique donors"),
        ("B", "Preparation relationship", "hDRG3 and hDRG5 were separate\npreparations from one unique donor"),
        ("C", "Locked detection gate", "Assessable: ≥20 nuclei per\npreparation/class\nPositive: fraction detected ≥0.05\nDonor positive if any assessable\npreparation was positive"),
    ]
    for j, (letter, title, txt) in enumerate(texts):
        ax = fig.add_subplot(gs[0, j]); panel(ax, letter, title); ax.axis("off"); ax.text(.06, .76, txt, va="top", fontsize=6.7, linespacing=1.35)

    ax = fig.add_subplot(gs[1, 0]); panel(ax, "D", "DRG candidate expression")
    z = expr[expr.gene.isin(genes) & expr.cell_type.str.lower().isin([c.lower() for c in classes])].copy()
    gx = {g: i for i, g in enumerate(classes)}; gy = {g: i for i, g in enumerate(genes)}
    for _, r in z.iterrows():
        canonical = {c.lower(): c for c in classes}[str(r.cell_type).lower()]
        ax.scatter(gx[canonical], gy[r.gene], s=max(3, float(r.fraction_expressing) * 260), c=float(r.mean_log_expression), cmap="viridis", vmin=z.mean_log_expression.min(), vmax=z.mean_log_expression.max())
    ax.set_xticks(range(len(classes)), ["Nociceptor", "Mechanoreceptor", "Schwann"], rotation=45, ha="right")
    ax.set_yticks(range(len(genes)), genes); ax.invert_yaxis(); tidy(ax, None)

    ax = fig.add_subplot(gs[1, 1]); panel(ax, "E", "Expression-positive recurrence")
    rr = rec[rec.gene.isin(genes) & rec.cell_type.str.lower().isin([c.lower() for c in classes])].copy()
    mat = np.full((len(genes), len(classes)), np.nan); labels = np.full(mat.shape, "", dtype=object)
    for _, r in rr.iterrows():
        canonical = {c.lower(): c for c in classes}[str(r.cell_type).lower()]
        i, j = genes.index(r.gene), classes.index(canonical)
        den = int(r.assessable_unique_donors); num = int(r.expression_positive_unique_donors)
        mat[i, j] = num / den if den else np.nan; labels[i, j] = f"{num}/{den}" if den else "NA"
    ax.imshow(np.ma.masked_invalid(mat), cmap="YlGnBu", vmin=0, vmax=1, aspect="auto")
    for i in range(len(genes)):
        for j in range(len(classes)):
            ax.text(j, i, labels[i, j], ha="center", va="center", fontsize=6.2, color="white" if np.isfinite(mat[i, j]) and mat[i, j] > .58 else DARK)
    ax.set_xticks(range(len(classes)), ["Nociceptor", "Mechanoreceptor", "Schwann"], rotation=45, ha="right")
    ax.set_yticks(range(len(genes)), genes)
    ax.text(.5, -.28, "positive / assessable unique donors", transform=ax.transAxes, ha="center", fontsize=5.9)

    ax = fig.add_subplot(gs[1, 2]); panel(ax, "F", "Threshold sensitivity in sensory neurons")
    ss = sens[sens.gene.isin(genes) & sens.cell_type.str.lower().isin(["nociceptor", "mechanoreceptor"]) & sens.minimum_nuclei.eq(20)].copy()
    for gene, g in ss.groupby("gene"):
        q = g.groupby("detection_fraction_threshold", as_index=False)["positive_unique_donors"].mean()
        ax.plot(q.detection_fraction_threshold, q.positive_unique_donors, marker="o", ms=3, lw=1, label=gene, color={"TXNL1":TX,"MAPK3":BLUE,"FGFR3":"#D8A62A","PDPR":"#817591","GFPT1":GREY}[gene])
    ax.axvline(.05, c=DARK, ls="--", lw=.7); ax.set(xlabel="Fraction-detected threshold", ylabel="Mean positive unique donors\nacross sensory classes")
    ax.legend(frameon=False, ncol=2, loc="best"); tidy(ax)
    expr.to_csv(SD / "S8_D_expression.tsv", sep="\t", index=False)
    rec.to_csv(SD / "S8_E_recurrence.tsv", sep="\t", index=False)
    sens.to_csv(SD / "S8_F_threshold_sensitivity.tsv", sep="\t", index=False)
    return save(fig, 8)


def s9():
    reg = load_fig(6, "G_H_registry")
    rec = pd.read_csv(PHASE / "03_drg_recurrence_fix/DRG_expression_positive_recurrence.tsv", sep="\t")
    genes = ["TXNL1", "MAPK3", "FGFR3", "PDPR", "GFPT1"]
    sensory = rec[rec.cell_type.str.lower().isin(["nociceptor", "mechanoreceptor"])].groupby("gene", as_index=False).agg(pos=("expression_positive_unique_donors", "max"), assess=("assessable_unique_donors", "max"))
    z = reg.set_index("gene").reindex(genes).reset_index().merge(sensory, on="gene", how="left")
    fig, aa = plt.subplots(2, 3, figsize=(7.0866, 5.4), gridspec_kw={"height_ratios":[1.15,.85]})
    ax = aa[0, 0]; panel(ax, "A", "Evidence domains")
    vals = np.column_stack([
        -np.log10(np.clip(pd.to_numeric(z.SCIATICA_MAGMA_FDR, errors="coerce"), 1e-12, 1)),
        np.abs(pd.to_numeric(z.robust_TWAS_z, errors="coerce")),
        pd.to_numeric(z.coloc_PP4, errors="coerce") * 6,
        pd.to_numeric(z.susie_coloc_max_PP4, errors="coerce") * 6,
        pd.to_numeric(z.gsMap_positive_sections, errors="coerce") / 20 * 6,
        z.pos / z.assess * 6,
    ])
    ax.imshow(np.ma.masked_invalid(vals), cmap="YlGnBu", aspect="auto")
    ax.set_yticks(range(len(genes)), genes); ax.set_xticks(range(6), ["MAGMA", "TWAS", "ABF", "SuSiE", "Spatial", "DRG+"], rotation=55, ha="right")
    for i, r in z.iterrows():
        ax.text(5, i, f"{int(r.pos)}/{int(r.assess)}" if pd.notna(r.assess) and r.assess else "NA", ha="center", va="center", fontsize=5.7, color="white" if pd.notna(r.pos) and r.pos/r.assess > .6 else DARK)
    ax = aa[0, 1]; panel(ax, "B", "Frozen candidate order")
    ax.axis("off"); ax.text(.05, .82, "1  TXNL1 — primary\n2  MAPK3 — backup\n3  FGFR3 — third", fontsize=10, color=SCI, linespacing=1.6, va="top")
    ax.text(.05, .25, "Corrected DRG recurrence did not\ncreate a hard conflict with the\nfrozen ranking.", fontsize=5.9, va="top", linespacing=1.3)
    ax = aa[0, 2]; panel(ax, "C", "Missing-evidence rule")
    ax.axis("off"); ax.text(.05, .82, "Unavailable evidence remained unavailable.\nSparse FGFR3 DRG detection was retained.\nNo composite numeric vote was used.", fontsize=7.2, linespacing=1.6, va="top")
    bottom = [
        ("D", "Leave-one-domain-out", "TXNL1 remained first when one\nevidence domain was omitted."),
        ("E", "Prespecified-weight sensitivity", "Alternative locked weights did not\nreverse the top candidate order."),
        ("F", "Claim boundary", "Prioritized candidate gene, not a\ncausal gene or therapeutic target."),
    ]
    for ax, (letter, title, txt) in zip(aa[1], bottom):
        panel(ax, letter, title); ax.axis("off"); ax.text(.04, .75, txt, fontsize=6.0, va="top", linespacing=1.35)
    z.to_csv(SD / "S9_candidate_evidence_corrected.tsv", sep="\t", index=False)
    return save(fig, 9)
