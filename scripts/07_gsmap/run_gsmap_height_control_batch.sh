#!/usr/bin/env bash
set -euo pipefail
root=.
cd "$root"
for sample in GSM6919905 GSM6919909 GSM6919911 GSM6919917; do
  scripts/common/phase3_attempt.sh gsmap "gsmap_height_control_${sample}_01" negative_control "GSE222322_${sample}" HEIGHT_IRN "LD_preserving_trait_control_QC_selected_one_per_donor" results/phase3_final/08_gsmap scripts/07_gsmap/run_gsmap_height_control_section.sh "$sample" &
done
wait
