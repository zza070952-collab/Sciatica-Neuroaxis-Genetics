#!/usr/bin/env python3
"""Prepare real, harmonized Phase 3 GWAS inputs for MAGMA/colocalization/gsMap."""
from pathlib import Path
import gzip
import pandas as pd

ROOT = Path(".")
IN = ROOT / "work/gwas_qc"
OUT = ROOT / "work/phase3_final/gwas_formats"
OUT.mkdir(parents=True, exist_ok=True)
traits = {
    "SCIATICA": (IN / "SCIATICA_R13.common_qc.tsv.gz", 28094, 347768),
    "LDH": (IN / "LDH_R13.common_qc.tsv.gz", 36573, 347768),
}
summary = []
for trait, (src, cases, controls) in traits.items():
    n_eff = 4.0 / (1.0/cases + 1.0/controls)
    magma = OUT / f"{trait}.magma.pval"
    magma_no_mhc = OUT / f"{trait}.magma.no_mhc.pval"
    meta = OUT / f"{trait}.coloc_variants.tsv.gz"
    gs = OUT / f"{trait}.gsmap.sumstats.gz"
    gs_no_mhc = OUT / f"{trait}.gsmap.no_mhc.sumstats.gz"
    first = True
    n = n_no_mhc = 0
    for c in pd.read_csv(src, sep="\t", chunksize=500_000):
        c = c.loc[c.chr.between(1, 22)].copy()
        if c.empty:
            continue
        c[["rsid", "p"]].rename(columns={"rsid": "SNP", "p": "P"}).to_csv(
            magma, sep="\t", index=False, mode="w" if first else "a", header=first)
        c["varID"] = "chr" + c.chr.astype(str) + "_" + c.position.astype(str) + "_" + c.ref + "_" + c.effect_allele + "_b38"
        c[["varID", "rsid", "chr", "position", "ref", "effect_allele", "beta", "se", "p"]].to_csv(
            meta, sep="\t", index=False, mode="wt" if first else "at", header=first,
            compression="gzip")
        g = pd.DataFrame({"SNP": c.rsid, "A1": c.effect_allele, "A2": c.ref,
                          "Z": c.z, "N": cases + controls})
        g.to_csv(gs, sep="\t", index=False, mode="wt" if first else "at", header=first,
                 compression="gzip")
        keep = ~((c.chr == 6) & c.position.between(25_000_000, 34_000_000))
        c.loc[keep, ["rsid", "p"]].rename(columns={"rsid": "SNP", "p": "P"}).to_csv(
            magma_no_mhc, sep="\t", index=False, mode="w" if first else "a", header=first)
        g.loc[keep].to_csv(gs_no_mhc, sep="\t", index=False,
                          mode="wt" if first else "at", header=first, compression="gzip")
        n += len(c); n_no_mhc += int(keep.sum()); first = False
    summary.append({"trait": trait, "autosomal_rows": n, "no_mhc_rows": n_no_mhc,
                    "n_case": cases, "n_control": controls, "effective_n": n_eff})
pd.DataFrame(summary).to_csv(OUT / "format_summary.tsv", sep="\t", index=False)
print(pd.DataFrame(summary).to_string(index=False))
