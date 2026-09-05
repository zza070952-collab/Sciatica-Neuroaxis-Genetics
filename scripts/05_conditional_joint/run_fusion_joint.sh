#!/usr/bin/env bash
set -euo pipefail
project=.
cd "$project"
thresholds=results/phase3_final/05_twas/fusion_bonferroni_thresholds.tsv
[[ -s "$thresholds" ]] || { echo "missing thresholds: $thresholds" >&2; exit 1; }
jobs=$(mktemp logs/phase3_final/twas/fusion_joint_jobs.XXXXXX)
trap 'rm -f "$jobs"' EXIT
for trait in SCIATICA LDH; do
  z=$(awk -F '\t' -v trait="$trait" 'NR>1 && $1==trait {print $4}' "$thresholds")
  [[ -n "$z" ]] || { echo "missing z threshold for $trait" >&2; exit 1; }
  for chr in $(seq 1 22); do printf '%s\t%s\t%s\n' "$trait" "$chr" "$z" >> "$jobs"; done
done
run_joint_logged() {
  trait=$1; chr=$2; z=$3
  scripts/common/phase3_attempt.sh twas "fusion_joint_${trait}_chr${chr}" primary "$trait" FUSION_5tissue_joint \
    "chr${chr};minp0.01;global_bonferroni_z=${z}" "results/phase3_final/05_twas/fusion_joint/${trait}" \
    scripts/05_conditional_joint/run_fusion_joint_one.sh "$trait" "$chr" "$z"
}
export -f run_joint_logged
export project
xargs -P 10 -n 3 bash -c 'cd "$project"; run_joint_logged "$1" "$2" "$3"' _ < "$jobs"
