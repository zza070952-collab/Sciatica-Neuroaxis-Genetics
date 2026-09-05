# Methods to code

| Final computational Methods | Production code | Coverage |
|---|---|---|
| GWAS acquisition/harmonization | scripts/01_gwas_qc/audit_and_prepare_gwas.py; phase3_gwas_qc.py; prepare_phase3_gwas_formats.py | common QC complete; production HM3 inputs inherited |
| LD score regression | scripts/02_ldsc/run_phase3_ldsc.sh; parse_phase3_ldsc.py | observed/liability h2 and rg |
| MAGMA genes/pathways | scripts/03_magma/run_phase3_magma.sh; run_magma_noMHC_pathways.sh; parse_phase3_magma.py | 10/10 kb, effective N, global correction, no-MHC |
| Five-tissue FUSION | scripts/04_fusion_twas/run_fusion_full.sh; run_fusion_one.sh; parse_phase3_fusion.py | two traits, five tissues, chr1-22 |
| Conditional/joint | scripts/05_conditional_joint/run_fusion_joint.sh; run_fusion_joint_one.sh; summarize_fusion_twas.py | production threshold and joint outputs |
| Pre-spatial genetics-only selection | scripts/06_coloc_susie/select_phase3_genetic_candidates.py | lookup dependencies explicitly declared |
| Full-cis ABF | scripts/06_coloc_susie/fetch_phase3_coloc_inputs.py; run_phase3_coloc_abf.R | seven QTD contexts, matched cis SNPs, priors |
| SuSiE colocalization | scripts/06_coloc_susie/prepare_phase3_susie_ld.py; run_phase3_susie_coloc.R | original EUR LD and allelic orientation |
| Adult spatial analysis | scripts/07_gsmap/prepare_phase3_gsmap_inputs.py; add_gsmap_pca_latent.py; run_gsmap_batch.sh; parse_phase3_gsmap.py | section/donor/ACAT/Fisher/LODO; corrected manifest required |
| Spatial sensitivities | scripts/07_gsmap/run_gsmap_no_mhc_sensitivity.sh; validate_gsmap_mhc_sensitivity.py; height-control scripts | retained MHC and reference-trait branches |
| Spinal/DRG references | scripts/08_spinal_snRNA/run_phase3_singlecell_base.py; run_phase3_spatial_label_transfer.py | inherited spinal annotation audit and DRG marker-based pipeline |
| DRG recurrence correction | scripts/09_DRG_snRNA/drg_recurrence.R | 6 preparations/5 donors; >=20 nuclei and >=5% detection |
| Integration and fixed nomination | scripts/10_candidate_integration/build_phase3_provisional_candidates.py; integrate_phase3_candidates_cells_space.py; freeze_phase3_candidates.py | fixed multi-domain order; no numerical vote |

The explicit upstream gaps are listed in `reproducibility_limitations.md`. Mapping a Methods heading to a retained consumer does not establish that its inherited inputs can be recreated from raw data by this release alone.
