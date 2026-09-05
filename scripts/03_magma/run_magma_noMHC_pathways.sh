#!/usr/bin/env bash
set -euo pipefail
root=.
cd "$root"
out=results/phase3_final/03_magma; magma=$(find work/software/phase3_magma/magma_v1.10 -type f -name magma -print -quit)
for trait in SCIATICA LDH; do
  "$magma" --gene-results "$out/raw/${trait}_10kb_noMHC.genes.raw" --set-annot raw/phase3_final/msigdb/c2.cp.v2026.1.Hs.entrez.gmt --out "$out/raw/${trait}_C2_noMHC"
  "$magma" --gene-results "$out/raw/${trait}_10kb_noMHC.genes.raw" --set-annot raw/phase3_final/msigdb/c5.go.v2026.1.Hs.entrez.gmt --out "$out/raw/${trait}_GO_noMHC"
done
