#!/usr/bin/env bash
set -euo pipefail

project=.
cd "$project"

manifest=work/software/phase3_twas/FUSION_WEIGHTS_EUR/archive_manifest.tsv
[[ -s "$manifest" ]] || { echo "missing completed archive manifest: $manifest" >&2; exit 1; }
[[ "$(awk 'NR>1 && $7=="PASS" {n++} END{print n+0}' "$manifest")" -eq 5 ]] || {
  echo 'not all five FUSION archives passed validation' >&2
  exit 1
}
if find work/software/phase3_twas/FUSION_WEIGHTS_EUR -name '*.aria2' -print -quit | grep -q .; then
  echo 'found aria2 sidecar in validated FUSION weights' >&2
  exit 1
fi

jobs=$(mktemp logs/phase3_final/twas/fusion_jobs.XXXXXX)
trap 'rm -f "$jobs"' EXIT
for trait in SCIATICA LDH; do
  for tissue in Nerve_Tibial Brain_Spinal_cord_cervical_c-1 Whole_Blood Muscle_Skeletal Cells_Cultured_fibroblasts; do
    for chr in $(seq 1 22); do
      printf '%s\t%s\t%s\n' "$trait" "$tissue" "$chr" >> "$jobs"
    done
  done
done

run_logged() {
  trait=$1
  tissue=$2
  chr=$3
  attempt="fusion_${trait}_${tissue}_chr${chr}"
  outdir="results/phase3_final/05_twas/fusion_raw/${trait}/${tissue}"
  outfile="$outdir/chr${chr}.dat"
  if [[ -s "$outfile" ]] && head -1 "$outfile" | grep -q $'^PANEL\tFILE\tID\tCHR\t'; then
    return 0
  fi
  scripts/common/phase3_attempt.sh twas "$attempt" primary "$trait" "GTExv8_EUR_${tissue}" \
    "chr${chr};clean_HM3;1000G_EUR_LD" "$outdir" \
    scripts/04_fusion_twas/run_fusion_one.sh "$trait" "$tissue" "$chr"
}
export -f run_logged
export project
xargs -P 20 -n 3 bash -c 'cd "$project"; run_logged "$1" "$2" "$3"' _ < "$jobs"

missing=0
for trait in SCIATICA LDH; do
  for tissue in Nerve_Tibial Brain_Spinal_cord_cervical_c-1 Whole_Blood Muscle_Skeletal Cells_Cultured_fibroblasts; do
    for chr in $(seq 1 22); do
      file="results/phase3_final/05_twas/fusion_raw/${trait}/${tissue}/chr${chr}.dat"
      if [[ ! -s "$file" ]] || [[ "$(head -1 "$file")" != PANEL$'\t'FILE$'\t'ID$'\t'CHR$'\t'P0$'\t'P1$'\t'HSQ$'\t'BEST.GWAS.ID$'\t'BEST.GWAS.Z$'\t'EQTL.ID$'\t'EQTL.R2$'\t'EQTL.Z$'\t'EQTL.GWAS.Z$'\t'NSNP$'\t'NWGT$'\t'MODEL$'\t'MODELCV.R2$'\t'MODELCV.PV$'\t'TWAS.Z$'\t'TWAS.P ]]; then
        printf 'invalid product: %s\n' "$file" >&2
        missing=$((missing+1))
      fi
    done
  done
done
[[ "$missing" -eq 0 ]]
