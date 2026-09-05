# Data acquisition and placement

Do not commit input data. Official acquisition URLs are in `../manifests/data_manifest.tsv`.

- FinnGen endpoint files: `data/gwas/FinnGen_DF13_M13_SCIATICA.gz` and `data/gwas/FinnGen_DF13_M13_LUMBAR_PROLAPSE.gz`. Obtain with author-approved FinnGen access. Height control, if used: `data/gwas/finngen_R13_HEIGHT_IRN.gz`.
- Production HapMap3 summaries: `work/gwas/M13_SCIATICA.sumstats.gz` and `work/gwas/M13_LUMBAR_PROLAPSE.sumstats.gz`; provide exact audited versions.
- LDSC EUR scores/weights: `raw/ld/ldsc_eur_w_ld_chr/` (chr1-22); LDSC source: `work/software/src/ldsc/`.
- MAGMA binary: `work/software/phase3_magma/magma_v1.10/`; EUR bed/bim/fam: `work/software/phase3_magma/g1000_eur/`; NCBI37.3 loc file: `work/software/phase3_magma/NCBI37.3/`.
- Locked MSigDB: `raw/phase3_final/msigdb/c2.cp.v2026.1.Hs.entrez.gmt` and `c5.go.v2026.1.Hs.entrez.gmt`; follow the resource's license/access terms.
- FUSION source: `work/software/phase3_twas/fusion_twas/`; LD prefix: `work/software/phase3_twas/LDREF/1000G.EUR.`; weights: `work/software/phase3_twas/FUSION_WEIGHTS_EUR/`, retaining the validated archive manifest. The acquisition-only helper is `scripts/04_fusion_twas/download_fusion_weights.sh`.
- GTEx PredictDB lookup databases: `work/software/phase3_twas/predictdb_core/eqtl/mashr/`; these supply gene identifiers/coordinates, not an additional association method.
- eQTL Catalogue full-cis data: scripts use tabix-accessible official URLs and write only queried locus data under `work/phase3_final/coloc/inputs/`.
- GSE222322: `raw/geo/GSE222322/unpacked/raw/`; corrected section metadata: `results/spatial_audit/spatial_manifest.tsv`.
- gsMap resources: `work/phase3_final/gsmap_resource/`, including the GENCODE lift37 gene annotation needed for gene lookup.
- GSE243076 prepared reference: `work/inherited/GSE243076_spinal_snRNA_results5s_annotated.h5ad`. See limitations for upstream-annotation provenance requirements.
- GSE168243: `raw/geo/GSE168243/unpacked/GSM*_prep*.csv.gz`; corrected computational S11B: `data/processed/S11B_DRG_expression.tsv`.
- Author computational Source Data: `data/computational_source_data/Source_Data_Figure1.xlsx` through `Source_Data_Figure6.xlsx`, using original sheet names. Full computational TSVs: `data/computational_source_tables/`.

Third-party matrices, weights, LD references and raw GWAS are not redistributed here. Accession access is not a substitute for preserving the exact release and checksums.
