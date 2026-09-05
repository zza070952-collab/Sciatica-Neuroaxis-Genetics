# Input contracts

All paths are relative to the checkout. Large inputs are excluded from Git.

| Stage | Inputs | Outputs |
|---|---|---|
| QC | FinnGen DF13 gzip TSV: #chrom, pos, ref, alt, rsids, pval, beta, sebeta, af_alt; original HapMap3 summaries | work/gwas_qc/*common_qc.tsv.gz, audit.json, Phase 3 QC tables and formats |
| LDSC | two production HapMap3 summaries SNP/A1/A2/Z/N; EUR LD scores and weights chr1-22 | observed/liability h2, intercept, ratio, rg |
| MAGMA | common-QC p-value files; MAGMA 1.10; NCBI37.3 gene locations; g1000 EUR; locked C2/GO GMT | primary/no-MHC gene and pathway results |
| FUSION | five validated GTExv8 EUR archives and pos/model files; 1000G EUR per-chromosome LD; harmonized SNP/A1/A2/Z/N | chr1-22 model results and global FDR |
| Joint | combined five-tissue FUSION .dat and same LD | retained conditional/joint models |
| Coloc | pre-spatial genetics-only candidates, exact full-cis QTD records and complete overlapping GWAS; PredictDB/GENCODE gene lookup | harmonized regions, PP0-PP4, SuSiE LD manifest and signal-pair results |
| gsMap | official section coordinate/count/image/scalefactor files, corrected donor manifest, gsMap 1.73.7 reference | section results, donor ACAT, across-donor ACAT/Fisher, LODO |
| Single-nucleus | annotated inherited GSE243076 H5AD and processed GSE168243 matrices | marker/QC/annotation audit and expression summaries |
| DRG correction | S11B TSV with gene, cell_type, donor, sample_preparation, biological_replicate_key, n_cells, fraction_expressing, mean_log_expression | assessable/positive preparation and unique-donor recurrence |
| Integration | frozen pre-spatial candidate table, core genetic/QTL results, spatial/cell evidence | preserved TXNL1/MAPK3/FGFR3 ranking |
| Figures | original computational Source Data workbooks; complete computational supplementary TSVs; corrected DRG tables | computational PDFs/SVGs/TIFFs and per-panel TSVs |

FUSION tissues: Nerve_Tibial, Brain_Spinal_cord_cervical_c-1, Whole_Blood, Muscle_Skeletal, Cells_Cultured_fibroblasts. Case and spelling match production filenames.

Main figure workbook sheets are defined exactly by `SHEET_MAP` in `scripts/11_figures_computational/main_figures.py`; Figure 5 also needs B_E_maps. The supplied author workbooks use A_FUSION/C_top_multisnp/D_joint/E_TXNL1_locus/F_coloc/G_context_coverage/B_candidate_heatmap and the corrected Figure 5/6 panel labels. Final exported panel TSVs have a different sheet naming scheme; do not silently substitute them.

S1-S4/S7 read the original computational TSV views under `data/computational_source_tables/Figure2/`, `Figure3/`, `Figure4/`, `Figure6/` and `07_supplementary_tables/`. S6/S8/S9 read the computational workbooks and corrected recurrence outputs. S5 uses the complete `spot_level_results.tsv.gz` and section summaries, not representative-section screenshots.

Source Data export is built into the plotting scripts (`write_source`/`write_source(n,...)` and the supplementary functions). The export does not require laboratory data.
