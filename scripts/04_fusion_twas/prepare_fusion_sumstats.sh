#!/usr/bin/env bash
set -euo pipefail

project=.
outdir="$project/work/phase3_final/gwas_formats"
mkdir -p "$outdir"
summary="$outdir/fusion_sumstats_qc.tsv"
printf 'trait\tinput_rows\tcomplete_rows\tdropped_incomplete\n' > "$summary"

for spec in SCIATICA:M13_SCIATICA LDH:M13_LUMBAR_PROLAPSE; do
  trait="${spec%%:*}"
  stem="${spec##*:}"
  input="$project/work/gwas/${stem}.sumstats.gz"
  output="$outdir/${trait}.fusion.sumstats.gz"
  input_rows="$(zcat "$input" | awk 'END {print NR-1}')"
  zcat "$input" | awk 'BEGIN {FS=OFS="\t"}
    NR==1 {print "SNP","A1","A2","Z","N"; next}
    NF==5 && $1!="" && $2~/^[ACGT]$/ && $3~/^[ACGT]$/ && $4~/^-?[0-9]+([.][0-9]+)?([eE][-+]?[0-9]+)?$/ && $5~/^[0-9]+([.][0-9]+)?$/ {print $1,$2,$3,$4,$5}' \
    | gzip -c > "$output"
  complete_rows="$(zcat "$output" | awk 'END {print NR-1}')"
  printf '%s\t%s\t%s\t%s\n' "$trait" "$input_rows" "$complete_rows" "$((input_rows-complete_rows))" >> "$summary"
  md5sum "$output" > "$output.md5"
  sha256sum "$output" > "$output.sha256"
done

md5sum "$summary" > "$summary.md5"
sha256sum "$summary" > "$summary.sha256"
