#!/usr/bin/env bash
set -euo pipefail
root=.
cd "$root"
sample=GSM6919905
# Retained successful production path: independent quick_mode, not the
# superseded run_spatial_ldsc call that had no retained annotation chunks.
bash scripts/07_gsmap/run_gsmap_section.sh "$sample" no_mhc
python scripts/07_gsmap/validate_gsmap_mhc_sensitivity.py
