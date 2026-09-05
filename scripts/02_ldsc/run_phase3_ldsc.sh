#!/usr/bin/env bash
set -euo pipefail

root=.
cd "$root"
out=results/phase3_final/02_genetic_architecture
mkdir -p "$out/ldsc_logs"
py=${LDSC_PYTHON:-python}
ldsc=work/software/src/ldsc/ldsc.py
ref=raw/ld/ldsc_eur_w_ld_chr/

run_h2() {
  local trait=$1 input=$2
  "$py" "$ldsc" --h2 "$input" --ref-ld-chr "$ref" --w-ld-chr "$ref" \
    --out "$out/ldsc_logs/${trait}_observed"
}

run_liability() {
  local trait=$1 input=$2 samp=$3 pop=$4 label=$5
  "$py" "$ldsc" --h2 "$input" --ref-ld-chr "$ref" --w-ld-chr "$ref" \
    --samp-prev "$samp" --pop-prev "$pop" --out "$out/ldsc_logs/${trait}_liability_${label}"
}

run_h2 SCIATICA work/gwas/M13_SCIATICA.sumstats.gz
run_h2 LDH work/gwas/M13_LUMBAR_PROLAPSE.sumstats.gz

# Population prevalences are sensitivity assumptions, not asserted epidemiologic truths.
run_liability SCIATICA work/gwas/M13_SCIATICA.sumstats.gz 0.074755 0.05 K005
run_liability SCIATICA work/gwas/M13_SCIATICA.sumstats.gz 0.074755 0.10 K010
run_liability LDH work/gwas/M13_LUMBAR_PROLAPSE.sumstats.gz 0.095160 0.03 K003
run_liability LDH work/gwas/M13_LUMBAR_PROLAPSE.sumstats.gz 0.095160 0.05 K005

"$py" "$ldsc" --rg work/gwas/M13_SCIATICA.sumstats.gz,work/gwas/M13_LUMBAR_PROLAPSE.sumstats.gz \
  --ref-ld-chr "$ref" --w-ld-chr "$ref" --out "$out/ldsc_logs/SCIATICA_LDH_rg"

python3 scripts/02_ldsc/parse_phase3_ldsc.py
