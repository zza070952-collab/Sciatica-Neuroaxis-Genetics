#!/usr/bin/env bash
set -euo pipefail
sample=$1
mode=${2:-main}
root=.
cd "$root"
config=scripts/07_gsmap/gsmap_sumstats_main.yaml
[[ "$mode" == no_mhc ]] && config=scripts/07_gsmap/gsmap_sumstats_no_mhc.yaml
gsmap quick_mode \
  --workdir "work/phase3_final/gsmap_runs/$mode" \
  --sample_name "$sample" \
  --gsMap_resource_dir work/phase3_final/gsmap_resource \
  --hdf5_path "work/phase3_final/gsmap_input_pca/${sample}.h5ad" \
  --annotation region --data_layer count \
  --latent_representation X_pca_phase3 \
  --sumstats_config_file "$config" --max_processes 8
