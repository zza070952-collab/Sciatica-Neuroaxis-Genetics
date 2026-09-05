#!/usr/bin/env bash
set -euo pipefail
sample=$1
root=.
cd "$root"
gsmap quick_mode \
  --workdir work/phase3_final/gsmap_runs/height_control --sample_name "$sample" \
  --gsMap_resource_dir work/phase3_final/gsmap_resource \
  --hdf5_path "work/phase3_final/gsmap_input_pca/${sample}.h5ad" \
  --annotation region --data_layer count --latent_representation X_pca_phase3 \
  --sumstats_config_file scripts/07_gsmap/gsmap_sumstats_height_control.yaml --max_processes 8
