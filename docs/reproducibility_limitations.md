# Exact scope and retained-input boundaries

This repository preserves the retained claim-bearing Phase 3 computational scripts and Phase 8 corrections. It does not claim a newly executed full-data rerun.

1. The original precomputed HapMap3 GWAS files enter LDSC/FUSION as explicit external prepared inputs. Their initial generating script was not recovered in the retained production package. The common-biallelic QC producer is included and is distinct from this HapMap3 input. Do not state that the included QC script alone reproduces the production HapMap3 summaries.
2. Spinal-cord analysis starts from the inherited annotated GSE243076 H5AD. The included scripts audit markers/QC and map this reference to Visium. The upstream raw-matrix integration/annotation recipe is not fully recovered; supplying arbitrary reprocessed annotations would not reproduce the locked analysis.
3. The spatial-manifest donor correction is an input to the section-preparation script. The manifest must match the official donor mapping used by the paper.
4. The final S5 PDF used an existing all-section vector atlas with crop/assembly adjustments. The included production `section_grid` renderer recovers all-section scientific content from the complete spot table; exact final pixel layout is not asserted.
5. S9's leave-one-domain-out/weight-sensitivity annotations are inherited text in the final figure script. No corresponding standalone numerical sensitivity producer was identified in the retained package. The release does not invent one or treat the annotations as a new computation.
6. Frozen R package versions are recorded where measured. The retained renv template lacked dependency version fields and is not represented as a complete lock file; `R_packages.tsv` is supplied instead. Environment YAMLs are specifications and were not all freshly solved during release preparation.
7. Figure 1's wet-lab card and all wet-lab analysis functions are omitted from this computational-only derivative. Full Figure 1 assembly remains in the author's private manuscript package.

These limitations must be resolved or explicitly accepted by the authors before claiming an end-to-end raw-data reproduction. They do not change the original statistical outputs.
