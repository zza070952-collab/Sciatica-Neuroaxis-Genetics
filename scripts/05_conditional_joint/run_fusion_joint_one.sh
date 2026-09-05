#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo 'usage: run_fusion_joint_one.sh TRAIT CHR ZTHRESH' >&2
  exit 2
fi
trait=$1
chr=$2
zthresh=$3
project="$(pwd)"
outdir="$project/results/phase3_final/05_twas/fusion_joint/$trait"
mkdir -p "$outdir"
cd "$project/work/software/phase3_twas/fusion_twas"
Rscript FUSION.post_process.R \
  --input "../../../../results/phase3_final/05_twas/fusion_combined_input/${trait}.dat" \
  --sumstats "../../../phase3_final/gwas_formats/${trait}.fusion.sumstats.gz" \
  --ref_ld_chr ../LDREF/1000G.EUR. \
  --chr "$chr" \
  --minp_input 0.01 \
  --zthresh "$zthresh" \
  --out "../../../../results/phase3_final/05_twas/fusion_joint/${trait}/chr${chr}"
