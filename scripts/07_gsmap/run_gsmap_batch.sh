#!/usr/bin/env bash
set -euo pipefail
root=.
cd "$root"
max_parallel=${1:-4}
samples=$(cut -f1 results/phase3_final/08_gsmap/section_qc.tsv | tail -n +2)
running=0
for sample in $samples; do
  scripts/common/phase3_attempt.sh gsmap "gsmap_main_${sample}" primary "GSE222322_${sample}" \
    gsMap_quick_mode PCA50_SCIATICA_LDH results/phase3_final/08_gsmap \
    scripts/07_gsmap/run_gsmap_section.sh "$sample" main &
  running=$((running+1))
  if [[ $running -ge $max_parallel ]]; then
    wait -n || true
    running=$((running-1))
  fi
done
wait
