#!/usr/bin/env python3
"""Create one independent AnnData object per GSE222322 Visium section.

No sections or donors are merged.  Coordinates are the official Space Ranger
tables previously gated in Phase 0/1.  The count layer remains raw integer
counts for gsMap.
"""
from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

ROOT = Path(".")
RAW = ROOT / "raw/geo/GSE222322/unpacked/raw"
MANIFEST = ROOT / "results/spatial_audit/spatial_manifest.tsv"
OUT = ROOT / "results/phase3_final/08_gsmap"
WORK = ROOT / "work/phase3_final/gsmap_input"


def digest(path: Path, kind: str) -> str:
    h = hashlib.new(kind)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def locate(gsm: str, suffix: str) -> Path:
    matches = list(RAW.glob(f"{gsm}_*_{suffix}"))
    if len(matches) != 1:
        raise RuntimeError(f"expected one {suffix} for {gsm}, found {matches}")
    return matches[0]


OUT.mkdir(parents=True, exist_ok=True)
WORK.mkdir(parents=True, exist_ok=True)
manifest = pd.read_csv(MANIFEST, sep="\t")
manifest = manifest.loc[manifest.dataset.eq("GSE222322")].copy()
if len(manifest) != 20 or manifest.gate.ne("PASS").any():
    raise RuntimeError("GSE222322 must contain exactly 20 coordinate-PASS sections")

rows = []
hash_rows = []
for rec in manifest.itertuples(index=False):
    h5 = locate(rec.sample, "filtered_feature_bc_matrix.h5")
    pos = locate(rec.sample, "tissue_positions.csv.gz")
    section_code = h5.name.removeprefix(rec.sample + "_").removesuffix("_filtered_feature_bc_matrix.h5")
    sf = ROOT / f"raw/geo/GSE222322/unpacked/scalefactors/scalefactors/{section_code}_scalefactors_json.json"
    image = Path(rec.image)
    if not sf.exists() or not image.exists():
        raise RuntimeError(f"missing scalefactor/image for {rec.sample}: {sf}, {image}")

    a = sc.read_10x_h5(h5)
    a.var_names_make_unique()
    a.obs_names = a.obs_names.astype(str)
    p = pd.read_csv(
        pos, header=None,
        names=["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"],
    ).set_index("barcode")
    common = a.obs_names.intersection(p.index)
    if len(common) != a.n_obs:
        raise RuntimeError(f"{rec.sample}: coordinate match {len(common)}/{a.n_obs}")
    a = a[common].copy()
    p = p.loc[a.obs_names]
    for col in p.columns:
        a.obs[col] = p[col].to_numpy()
    a.obs["sample"] = rec.sample
    a.obs["section"] = rec.section
    a.obs["donor"] = rec.donor
    a.obs["patient_id"] = str(rec.patient_id)
    a.obs["region"] = "adult lumbar spinal cord"
    a.obsm["spatial"] = p[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy(dtype=float)
    if not sparse.issparse(a.X):
        a.X = sparse.csr_matrix(a.X)
    a.layers["count"] = a.X.copy()
    a.var["mt"] = a.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], inplace=True, log1p=False, percent_top=None)
    a.uns["spatial_provenance"] = {
        "coordinate_file": str(pos), "scalefactors_file": str(sf),
        "image_file": str(image), "coordinate_gate": "PASS",
        "donor_mapping": "results/spatial_audit/GSE222322_corrected_donor_mapping.tsv",
    }
    with sf.open() as handle:
        a.uns["scalefactors"] = json.load(handle)
    out = WORK / f"{rec.sample}.h5ad"
    a.write_h5ad(out, compression="gzip")
    rows.append({
        "sample": rec.sample, "donor": rec.donor, "patient_id": rec.patient_id,
        "section": rec.section, "spots": a.n_obs, "genes": a.n_vars,
        "median_umi": float(np.median(a.obs.total_counts)),
        "median_genes": float(np.median(a.obs.n_genes_by_counts)),
        "median_pct_mt": float(np.median(a.obs.pct_counts_mt)),
        "coordinate_match_rate": len(common) / a.n_obs,
        "coordinate_gate": "PASS", "biological_unit": "donor",
        "h5ad": str(out),
    })
    hash_rows.append({"file": str(out.relative_to(ROOT)), "size": out.stat().st_size,
                      "md5": digest(out, "md5"), "sha256": digest(out, "sha256")})

pd.DataFrame(rows).to_csv(OUT / "section_qc.tsv", sep="\t", index=False)
pd.DataFrame(hash_rows).to_csv(OUT / "input_h5ad_checksums.tsv", sep="\t", index=False)
print(json.dumps({"sections": len(rows), "donors": len(set(x["donor"] for x in rows)),
                  "spots": sum(x["spots"] for x in rows)}, indent=2))
