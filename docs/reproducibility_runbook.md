# Reproducibility runbook

1. Verify `checksums.sha256` against this checkout before changes.
2. Read the limitations and the input contracts. Prepare resources at the declared relative locations, including the author-supplied inherited inputs and their source checksums.
3. Create separate environments; install upstream FUSION at the recorded commit and the recorded MAGMA/gsMap versions. Preserve upstream license terms. Use the production LD build and alleles.
4. Run `python workflow/run.py --dry-run` to identify missing files; it does not download data or execute analyses.
5. Run `bash tests/smoke_test.sh` with Python3, NumPy, PyYAML and R available.
6. Run stages in workflow manifest order. Inspect command logs and output tables before advancing. Failed commands raise nonzero exit status; no success flag is fabricated.
7. Apply the corrected DRG recurrence after the preparation/donor mapping is verified. The earlier raw donor-column values are preparation identifiers and must not be interpreted as six independent people.
8. Rebuild computational figures only from source-locked data. Compare source TSV numeric values rather than PDF byte hashes when system fonts/renderer versions differ.
9. Record environment exports, input checksums and run logs for any independent rerun. No complete rerun is claimed by the publication-preparation smoke test.
