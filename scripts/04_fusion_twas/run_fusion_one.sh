#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo 'usage: run_fusion_one.sh TRAIT TISSUE CHR' >&2
  exit 2
fi

trait=$1
tissue=$2
chr=$3
project="$(pwd)"
outdir="$project/results/phase3_final/05_twas/fusion_raw/$trait/$tissue"
mkdir -p "$outdir"

cd "$project/work/software/phase3_twas/fusion_twas"
Rscript FUSION.assoc_test.R \
  --sumstats "../../../phase3_final/gwas_formats/${trait}.fusion.sumstats.gz" \
  --weights "../FUSION_WEIGHTS_EUR/extracted/${tissue}/GTExv8.EUR.${tissue}.pos" \
  --weights_dir "../FUSION_WEIGHTS_EUR/extracted/${tissue}" \
  --ref_ld_chr ../LDREF/1000G.EUR. \
  --chr "$chr" \
  --out "../../../../results/phase3_final/05_twas/fusion_raw/${trait}/${tissue}/chr${chr}.dat"
