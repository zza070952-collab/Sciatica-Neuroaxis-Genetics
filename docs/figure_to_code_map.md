# Computational figures to code

| Figure | Renderer | Locked source |
|---|---|---|
| Figure 1, computational blocks | main_figures.py::fig1 | numerical foundation plus Figures 2/5/6 Source Data; wet-lab card excluded |
| Figure 2 | main_figures.py::fig2 | GWAS, loci, QQ, LDSC and rg workbook sheets |
| Figure 3 | main_figures.py::fig3 | MAGMA/pathway/reference-comparison workbook sheets |
| Figure 4 | main_figures.py::fig4 | FUSION/joint/full-cis/SuSiE workbook sheets |
| Figure 5 | main_figures.py::fig5 | section maps, donor/meta/LODO/reference-control/overlap sheets |
| Figure 6 | main_figures.py::fig6 | spinal/DRG/reference association and Phase 8 corrected recurrence |
| S1-S4 | supplementary_base.py::s1/s2/s3/s4 | complete computational Source Data TSVs |
| S5 | all_section_atlas.py::section_grid | Phase 3 all-spot and section results; final PDF crop not reproduced |
| S6 | supplementary_updates.py::s6 | donor/meta/LODO and sensitivity sheets |
| S7 | supplementary_base.py::s7_s8_s9(7) | spinal UMAP/composition/expression TSVs |
| S8-S9 | supplementary_updates.py::s8/s9 | original computational sheets plus Phase 8 corrected recurrence |

All figure scripts are in `scripts/11_figures_computational/`. `run_supplementary.py` provides the S1-S9 entry point. Each data renderer exports machine-readable Source Data; no low-resolution screenshot assembly is used. Exact full manuscript visual assembly is outside the computational-only scope. No laboratory results renderer is included.
