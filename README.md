# Sciatica Neuroaxis Genetics

Reproducible computational workflows supporting the manuscript on integrated genetic, spatial and sensory-neuroaxis prioritization of TXNL1 in sciatica.

This repository contains manuscript computational analyses and computational figure code. Wet-laboratory raw data, statistical scripts, and laboratory records are excluded. The primary phenotype is FinnGen DF13 M13_SCIATICA; M13_LUMBAR_PROLAPSE is a structural reference from the same release, not independent replication. The frozen candidate order is TXNL1, MAPK3, then FGFR3.

## Workflow and scope

GWAS QC → LDSC → MAGMA → five-tissue FUSION → conditional/joint → full-cis ABF/SuSiE colocalization → gsMap → adult spinal/DRG single-nucleus references → candidate integration → computational figures.

The source-to-release SHA256 map is in `manifests/production_source_map.tsv`. Production statistical logic is retained; paths, entry points, and mixed-workbook reads are adapted for release. The repository does not rerun excluded exploratory methods. This is a computational-source release, not a claim that a fresh full-data run has been executed during publication preparation.

## System requirements and installation

The production scripts target Linux x86-64, Bash, Python 3, R, Conda/Mamba, GNU core utilities and GNU time. Run commands from the repository root. Individual stages require separate environments. Production software versions and exact recorded commits are in `manifests/software_versions.tsv` and `R_packages.tsv`; missing version records are explicitly marked. The YAML files are installation specifications, not a fabricated complete historical lock file.

```bash
conda env create -f environment.yml
conda env create -f envs/qc.yaml
conda env create -f envs/ldsc.yaml
conda env create -f envs/magma.yaml
conda env create -f envs/twas.yaml
conda env create -f envs/coloc.yaml
conda env create -f envs/gsmap.yaml
conda env create -f envs/singlecell.yaml
conda env create -f envs/drg.yaml
conda env create -f envs/figures.yaml
```

Install the upstream executables and reference resources into the relative locations in `data/README.md`. LDSC 1.0.1 uses its separate legacy environment; downstream parsers need Python 3. The legacy script runner accepts `LDSC_PYTHON` for the LDSC executable. Do not resolve R/Python version conflicts by mixing stages into one environment.

## Data preparation

See `manifests/data_manifest.tsv`, `manifests/expected_inputs.tsv` and `docs/input_contracts.md` for resource URLs, releases, file names and schemas. Obtain FinnGen data through its official access process. GTEx/FUSION weights, LD panels, licensed pathway resources and GEO matrices are not redistributed. No bulk dataset is downloaded by the smoke test.

Three inherited inputs require particular care: the production HapMap3 summaries, the corrected Space Ranger section/donor manifest, and the annotated spinal-cord reference object. They were inputs to the frozen Phase 3 analysis. Their original upstream construction is not fully contained in the retained Phase 3 code. See `docs/reproducibility_limitations.md` before describing this archive as a de novo raw-data reproduction. These inputs must be supplied with the matching manifests/checksums; a missing input is an error, not an invitation to substitute an unverified file.

Main computational figure workbooks and complete computational supplementary TSVs are author Source Data, distributed separately with the manuscript. Place only computational workbooks under `data/computational_source_data/`. DRG recurrence takes a standalone S11B TSV and does not need the mixed supplementary workbook.

## Running

The unified entry point is `workflow/run.py`. It lists required inputs without executing scientific analyses:

```bash
conda run -n sciatica-workflow python workflow/run.py --dry-run
conda run -n sciatica-workflow python workflow/run.py --stage gwas_qc --conda
conda run -n sciatica-workflow python workflow/run.py --stage magma,fusion,joint --conda
conda run -n sciatica-workflow python workflow/run.py --stage coloc --conda
conda run -n sciatica-workflow python workflow/run.py --stage gsmap,singlecell,integration --conda
conda run -n sciatica-workflow python workflow/run.py --stage drg_recurrence,figures --conda
```

After preparing all input contracts, `--stage all --conda` runs the ordered stages. The stage manifest records commands and environments; scripts retain their original explicit relative reference locations. `config/config.example.yaml` documents frozen parameters rather than overriding unmodified constants in production code. `config/paths.example.yaml` documents the required layout.

Outputs include QC variant counts and harmonized summaries, LDSC heritability/correlation, MAGMA genes/pathways, all five-tissue FUSION and joint results, ABF/SuSiE tables, section/donor/meta/leave-one-donor-out gsMap results, single-nucleus expression/annotation, corrected DRG recurrence and the locked candidate registry. Main computational figures write PDF/SVG/TIFF plus panel Source Data TSVs. Mapping: `docs/methods_to_code_map.md` and `docs/figure_to_code_map.md`.

## Runtime and memory

The retained successful-production logs are summarized in `manifests/runtime_memory.tsv` when available. A missing measurement is `NOT_RECORDED`; no estimated number is asserted as measured. Full-run timing depends on data access, hardware and section concurrency. The release-preparation smoke test is not a genome-wide benchmark. FUSION uses the production 20-job concurrency and gsMap 8 processes per section; review concurrent memory needs before running on a smaller machine.

## Smoke test

```bash
conda activate sciatica-qc
bash tests/smoke_test.sh
```

R must also be available on PATH for this cross-language syntax test. This starts Python/R, parses configuration and scripts, resolves stage references, runs the actual production ACAT/BH functions on analytical software-test inputs, and writes `results/smoke/smoke_output.json`. It neither simulates research results nor downloads protected data. The recorded execution report is `docs/smoke_test_report.md`.

## Reproducibility and interpretation

The QC implementation reports the release's absent per-variant INFO and retains its original adjacent-duplicate and palindromic-SNP logic. LDH cannot nominate a sciatica candidate. Normal adult atlases are expression/localization references, not case-control tests. Spatial aggregation uses donors. GSE168243 comprises six preparations from five unique donors; use the Phase 8 recurrence correction for final recurrence counts. There are no validated dorsal/ventral labels in these spatial plots. Height is an LD-preserving calibration/reference trait, not a guaranteed null trait.

S9 contains inherited narrative annotations about ranking sensitivity; these should not be mistaken for a newly implemented numerical re-ranking analysis. No new ranking method is introduced by this release. Consult the limitations document for exact coverage boundaries.

## Troubleshooting

- Missing file: inspect the stage prerequisites and input manifest; obtain the matching source/version.
- FUSION: verify all five weight archives and per-chromosome LD before starting. Missing/partial archives fail the original gate.
- Colocalization: verify build, SNP IDs, alleles, complete cis coverage and LD dimensions. Missing QTL evidence stays unavailable.
- Environment mismatch: activate the stage-specific environment and compare recorded versions. Do not change thresholds to bypass numerical failures.
- Figure input schema: use the original panel sheet names listed in the input contracts, not a manually relabeled workbook.

## Citation and license

Repository: https://github.com/zza070952-collab/Sciatica-Neuroaxis-Genetics

See `CITATION.cff` and cite the original resource/software papers. No archival DOI has been minted. No author email is invented.

License decision pending: see `LICENSE_PENDING.md`. Public visibility by itself does not grant an MIT/Apache/GPL reuse license. A formal versioned release awaits the author's explicit license decision.
