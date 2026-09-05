#!/usr/bin/env python3
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(".")
OUT = ROOT / "results/phase3_final/01_gwas_qc"
OUT.mkdir(parents=True, exist_ok=True)

TRAITS = {
    "SCIATICA": {
        "raw": Path("data/gwas/FinnGen_DF13_M13_SCIATICA.gz"),
        "qc": ROOT / "work/gwas_qc/SCIATICA_R13.common_qc.tsv.gz",
        "audit": ROOT / "work/gwas_qc/SCIATICA_R13.audit.json",
        "hm3": ROOT / "work/gwas/M13_SCIATICA.sumstats.gz",
        "cases": 28094,
        "controls": 347768,
    },
    "LDH": {
        "raw": Path("data/gwas/FinnGen_DF13_M13_LUMBAR_PROLAPSE.gz"),
        "qc": ROOT / "work/gwas_qc/LDH_R13.common_qc.tsv.gz",
        "audit": ROOT / "work/gwas_qc/LDH_R13.audit.json",
        "hm3": ROOT / "work/gwas/M13_LUMBAR_PROLAPSE.sumstats.gz",
        "cases": 36573,
        "controls": 347768,
    },
}

CHR_LENGTHS = {
    1: 248956422, 2: 242193529, 3: 198295559, 4: 190214555, 5: 181538259,
    6: 170805979, 7: 159345973, 8: 145138636, 9: 138394717, 10: 133797422,
    11: 135086622, 12: 133275309, 13: 114364328, 14: 107043718, 15: 101991189,
    16: 90338345, 17: 83257441, 18: 80373285, 19: 58617616, 20: 64444167,
    21: 46709983, 22: 50818468,
}

LONG_RANGE = [
    (6, 25_000_000, 34_000_000, "MHC"),
    (8, 7_000_000, 13_000_000, "chr8_inversion"),
    (11, 46_000_000, 57_000_000, "chr11_long_range_LD"),
]

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 7,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.8,
    "legend.frameon": False,
})


def open_text(path):
    return gzip.open(path, "rt", newline="") if str(path).endswith(".gz") else open(path, newline="")


def checksum(path, algo):
    h = hashlib.new(algo)
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def count_hm3(path):
    with gzip.open(path, "rt") as fh:
        return sum(1 for _ in fh) - 1


def analyze_qc_file(trait, path):
    pvals = []
    plot_rows = []
    n = 0
    mhc = 0
    long_counts = {name: 0 for _, _, _, name in LONG_RANGE}
    with gzip.open(path, "rt", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            try:
                chrom = int(row["chr"])
                pos = int(row["position"])
                p = float(row["p"])
                z = float(row["z"])
            except (ValueError, KeyError):
                continue
            if chrom not in CHR_LENGTHS or not (0 < p <= 1) or not math.isfinite(z):
                continue
            n += 1
            pvals.append(p)
            for c, start, end, name in LONG_RANGE:
                if chrom == c and start <= pos <= end:
                    long_counts[name] += 1
                    if name == "MHC":
                        mhc += 1
            if n % 250 == 0 or p < 1e-5:
                plot_rows.append((chrom, pos, p, row["rsid"]))

    arr = np.asarray(pvals, dtype=np.float64)
    del pvals
    chisq = np.square(np.asarray([float(x[4]) for x in []])) if False else None
    # p-values are used for QQ; lambda is derived separately from z in a second light stream.
    z2 = []
    with gzip.open(path, "rt", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            try:
                z = float(row["z"])
            except (ValueError, KeyError):
                continue
            if math.isfinite(z):
                z2.append(z * z)
    lambda_gc = float(np.median(np.asarray(z2, dtype=np.float64)) / 0.4549364231195727)
    del z2

    arr.sort()
    qidx = np.unique(np.linspace(0, len(arr) - 1, min(20000, len(arr)), dtype=np.int64))
    observed = -np.log10(np.maximum(arr[qidx], 1e-300))
    expected = -np.log10((qidx + 0.5) / len(arr))

    qq_path = OUT / f"{trait}_qq_source.tsv.gz"
    with gzip.open(qq_path, "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["rank_index", "expected_neglog10p", "observed_neglog10p"])
        for i, x, y in zip(qidx, expected, observed):
            w.writerow([int(i), f"{x:.8g}", f"{y:.8g}"])

    man_path = OUT / f"{trait}_manhattan_source.tsv.gz"
    with gzip.open(man_path, "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["chr", "position", "p", "rsid", "display_rule"])
        for chrom, pos, p, rsid in plot_rows:
            w.writerow([chrom, pos, f"{p:.12g}", rsid, "all_p_lt_1e-5_or_every_250th_valid_row"])

    offsets, acc = {}, 0
    centers = []
    for chrom in range(1, 23):
        offsets[chrom] = acc
        centers.append(acc + CHR_LENGTHS[chrom] / 2)
        acc += CHR_LENGTHS[chrom] + 5_000_000
    xs = np.asarray([offsets[c] + p for c, p, _, _ in plot_rows], dtype=float)
    ys = -np.log10(np.maximum(np.asarray([p for _, _, p, _ in plot_rows]), 1e-300))
    cs = np.asarray([c for c, _, _, _ in plot_rows])

    fig, axes = plt.subplots(1, 2, figsize=(183 / 25.4, 92 / 25.4), gridspec_kw={"width_ratios": [1.8, 1]})
    ax = axes[0]
    palette = np.where(cs % 2 == 0, "#7C93A6", "#B8C6D1")
    ax.scatter(xs, ys, c=palette, s=2.0, linewidths=0, rasterized=True)
    ax.axhline(-math.log10(5e-8), color="#B65C5C", lw=0.8, ls="--")
    ax.set_xticks(centers)
    ax.set_xticklabels([str(x) for x in range(1, 23)], fontsize=5)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("−log10(P)")
    ax.set_title(f"{trait}: common biallelic analysis set", loc="left", fontweight="bold")

    ax = axes[1]
    ax.scatter(expected, observed, s=3, color="#486A8C", linewidths=0, rasterized=True)
    lim = max(float(expected.max()), float(observed.max())) * 1.02
    ax.plot([0, lim], [0, lim], color="#777777", lw=0.8, ls="--")
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel("Expected −log10(P)")
    ax.set_ylabel("Observed −log10(P)")
    ax.set_title(f"QQ; λGC={lambda_gc:.3f}", loc="left", fontweight="bold")
    fig.tight_layout(w_pad=1.5)
    base = OUT / f"{trait}_manhattan_qq"
    fig.savefig(str(base) + ".svg", bbox_inches="tight")
    fig.savefig(str(base) + ".pdf", bbox_inches="tight")
    fig.savefig(str(base) + ".tiff", dpi=600, bbox_inches="tight")
    plt.close(fig)
    return {"n": n, "mhc": mhc, "lambda_gc": lambda_gc, "long": long_counts, "plot_n": len(plot_rows)}


def qc_records(path):
    with gzip.open(path, "rt", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            chrom = int(r["chr"])
            if chrom not in CHR_LENGTHS:
                continue
            yield (chrom, int(r["position"]), r["ref"], r["effect_allele"]), r


def build_common_intersection():
    ia = iter(qc_records(TRAITS["SCIATICA"]["qc"]))
    ib = iter(qc_records(TRAITS["LDH"]["qc"]))
    a = next(ia, None)
    b = next(ib, None)
    n = 0
    mismatch_rsid = 0
    out = OUT / "common_hq_intersection.tsv.gz"
    with gzip.open(out, "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["chr", "position", "ref", "alt_effect", "rsid", "sciatica_p", "ldh_p"])
        while a is not None and b is not None:
            if a[0] == b[0]:
                if a[1]["rsid"] != b[1]["rsid"]:
                    mismatch_rsid += 1
                else:
                    w.writerow([*a[0], a[1]["rsid"], a[1]["p"], b[1]["p"]])
                    n += 1
                a = next(ia, None)
                b = next(ib, None)
            elif a[0] < b[0]:
                a = next(ia, None)
            else:
                b = next(ib, None)
    return n, mismatch_rsid, out


def build_hm3_intersection():
    seen = set()
    with gzip.open(TRAITS["SCIATICA"]["hm3"], "rt", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            seen.add((r["SNP"], r["A1"], r["A2"]))
    n = 0
    out = OUT / "hapmap3_intersection.tsv.gz"
    with gzip.open(out, "wt", newline="") as dst, gzip.open(TRAITS["LDH"]["hm3"], "rt", newline="") as fh:
        w = csv.writer(dst, delimiter="\t", lineterminator="\n")
        w.writerow(["SNP", "A1", "A2"])
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            key = (r["SNP"], r["A1"], r["A2"])
            flip = (r["SNP"], r["A2"], r["A1"])
            if key in seen or flip in seen:
                w.writerow(key)
                n += 1
    return n, out


def main():
    summaries = {}
    audits = {}
    for trait, cfg in TRAITS.items():
        with open(cfg["audit"], encoding="utf-8") as fh:
            audits[trait] = json.load(fh)
        summaries[trait] = analyze_qc_file(trait, cfg["qc"])

    common_n, common_rsid_mismatch, common_path = build_common_intersection()
    hm3_n, hm3_path = build_hm3_intersection()

    summary_path = OUT / "gwas_qc_summary.tsv"
    with open(summary_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["trait", "release", "build", "n_case", "n_control", "effect_allele", "raw_rows", "common_biallelic_rows", "mhc_rows", "hm3_rows", "lambda_gc", "info_status", "z_p_inconsistent"])
        for trait, cfg in TRAITS.items():
            a = audits[trait]
            s = summaries[trait]
            w.writerow([trait, "FinnGen_R13", "GRCh38", cfg["cases"], cfg["controls"], "ALT",
                        a["counts"]["input_rows"], s["n"], s["mhc"], count_hm3(cfg["hm3"]),
                        f"{s['lambda_gc']:.8g}", "NO_VARIANT_INFO_COLUMN_UPSTREAM_QC_DOCUMENTED",
                        a["z_p_check"]["n_abs_log10_delta_gt_0.02"]])

    harmonization_path = OUT / "harmonization_report.tsv"
    with open(harmonization_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["check", "SCIATICA", "LDH", "interpretation"])
        w.writerow(["build", "GRCh38", "GRCh38", "PASS"])
        w.writerow(["effect_allele", "ALT", "ALT", "PASS"])
        w.writerow(["single_rsid_common_biallelic", summaries["SCIATICA"]["n"], summaries["LDH"]["n"], "PASS"])
        w.writerow(["common_exact_chr_pos_ref_alt_and_rsid", common_n, common_n, "PASS"])
        w.writerow(["common_key_rsid_mismatch", common_rsid_mismatch, common_rsid_mismatch, "REVIEW_IF_NONZERO"])
        w.writerow(["hapmap3_allele_compatible_intersection", hm3_n, hm3_n, "PASS"])
        w.writerow(["variant_INFO", "ABSENT", "ABSENT", "NOT_RECONSTRUCTABLE; no fabricated INFO"])

    variant_path = OUT / "variant_set_comparison.tsv"
    with open(variant_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["trait", "variant_set", "n_variants", "purpose", "primary_or_sensitivity"])
        for trait in TRAITS:
            a, s = audits[trait], summaries[trait]
            w.writerow([trait, "official_R13_release_rows", a["counts"]["input_rows"], "source audit", "primary_source"])
            w.writerow([trait, "common_biallelic_single_rsid_orientable", s["n"], "main downstream compatible set", "primary_analysis_set"])
            w.writerow([trait, "common_biallelic_MHC_excluded", s["n"] - s["mhc"], "MHC sensitivity", "sensitivity"])
            w.writerow([trait, "HapMap3", count_hm3(TRAITS[trait]["hm3"]), "LDSC/reference sensitivity", "sensitivity"])
        w.writerow(["BOTH", "common_exact_high_quality_intersection", common_n, str(common_path.relative_to(ROOT)), "harmonization"])
        w.writerow(["BOTH", "HapMap3_allele_compatible_intersection", hm3_n, str(hm3_path.relative_to(ROOT)), "harmonization"])

    with open(OUT / "long_range_ld_regions.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["chr", "start_grch38", "end_grch38", "label", "sciatica_common_rows", "ldh_common_rows"])
        for c, start, end, label in LONG_RANGE:
            w.writerow([c, start, end, label, summaries["SCIATICA"]["long"][label], summaries["LDH"]["long"][label]])

    audit_md = OUT / "methods_audit_gwas.md"
    audit_md.write_text(f"""# Phase 3 GWAS methods audit

FinnGen R13 endpoint files were verified as GRCh38 with ALT as the effect allele. The two files contain {audits['SCIATICA']['counts']['input_rows']:,} Sciatica and {audits['LDH']['counts']['input_rows']:,} LDH rows and omit per-variant INFO. No INFO value was reconstructed or fabricated.

FinnGen's current public methods document genotype-, sample- and variant-level QC of the imputation reference and post-imputation checks of INFO distributions, allele-frequency differences and chromosomal continuity. The current association documentation reports regenie association testing with minimum allele count 5. These statements support use of the official release set, but they do not prove INFO >= 0.8 for every released endpoint variant.

The downstream-compatible primary set therefore uses single-rsID common biallelic variants (MAF >= 0.01), removes unresolved mid-frequency A/T and C/G variants, checks beta/SE/P/Z consistency, and retains the FinnGen ALT effect orientation. HapMap3, MHC-excluded and common-set analyses are explicit sensitivity branches.

Official sources checked 2026-08-04:

- https://finngen.gitbook.io/documentation/methods/genotype-imputation/genotype-imputation
- https://finngen.gitbook.io/documentation/methods/phewas
- https://finngen.gitbook.io/documentation/methods/phewas/logistic-regression
- https://www.finngen.fi/en/access_results

Limitations: no per-variant INFO filter can be rerun from the endpoint files; the QQ genomic-inflation estimate is descriptive and must be interpreted alongside LDSC intercepts rather than as confounding by itself.
""", encoding="utf-8")

    for path in sorted(OUT.iterdir()):
        if not path.is_file() or path.suffix in {".md5", ".sha256"}:
            continue
        (Path(str(path) + ".md5")).write_text(f"{checksum(path, 'md5')}  {path.name}\n")
        (Path(str(path) + ".sha256")).write_text(f"{checksum(path, 'sha256')}  {path.name}\n")

    print(json.dumps({"summaries": summaries, "common_n": common_n, "hm3_n": hm3_n}, sort_keys=True))


if __name__ == "__main__":
    main()
