#!/usr/bin/env bash
set -euo pipefail
root=.
cd "$root"
out=results/phase3_final/03_magma
mkdir -p "$out/raw" "$out/regional_plots"
magma=$(find work/software/phase3_magma/magma_v1.10 -type f -name magma -print -quit)
bim=$(find work/software/phase3_magma/g1000_eur -type f -name '*.bim' -print -quit)
base=${bim%.bim}
gene_loc=$(find work/software/phase3_magma/NCBI37.3 -type f -name '*.loc' -print -quit)
test -x "$magma"; test -s "$base.bed"; test -s "$gene_loc"

"$magma" --annotate window=10,10 --snp-loc "$bim" --gene-loc "$gene_loc" --out "$out/raw/primary_10kb"

run_trait() {
  local trait=$1 n=$2
  "$magma" --bfile "$base" --pval "work/phase3_final/gwas_formats/${trait}.magma.pval" N="$n" \
    --gene-annot "$out/raw/primary_10kb.genes.annot" --out "$out/raw/${trait}_10kb"
  "$magma" --gene-results "$out/raw/${trait}_10kb.genes.raw" \
    --set-annot raw/phase3_final/msigdb/c2.cp.v2026.1.Hs.entrez.gmt \
    --out "$out/raw/${trait}_C2"
  "$magma" --gene-results "$out/raw/${trait}_10kb.genes.raw" \
    --set-annot raw/phase3_final/msigdb/c5.go.v2026.1.Hs.entrez.gmt \
    --out "$out/raw/${trait}_GO"
  "$magma" --bfile "$base" --pval "work/phase3_final/gwas_formats/${trait}.magma.no_mhc.pval" N="$n" \
    --gene-annot "$out/raw/primary_10kb.genes.annot" --out "$out/raw/${trait}_10kb_noMHC"
}
# MAGMA requires integer sample sizes in the --pval N modifier.  Effective
# sample sizes are therefore rounded to the nearest integer for computation;
# the exact values remain recorded below for auditability.
run_trait SCIATICA 103976
run_trait LDH 132371

{
  echo -e 'parameter\tvalue'
  echo -e 'genome_build\tGRCh37 MAGMA reference; rsID harmonization from FinnGen GRCh38'
  echo -e 'gene_window\t10kb upstream and 10kb downstream'
  echo -e 'sciatica_effective_n\t103976.39763530232'
  echo -e 'ldh_effective_n\t132371.1918738828'
  echo -e 'magma_sciatica_N\t103976 (nearest integer)'
  echo -e 'magma_ldh_N\t132371 (nearest integer)'
  echo -e 'MHC_sensitivity\tchr6:25-34Mb SNPs excluded before gene analysis'
} > "$out/magma_parameters.tsv"
