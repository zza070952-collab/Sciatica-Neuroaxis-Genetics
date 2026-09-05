#!/usr/bin/env bash
set -euo pipefail
root=.
cd "$root"
mkdir -p work/phase3_final/coloc
src=work/phase3_final/gwas_formats/SCIATICA.coloc_variants.tsv.gz
dst=work/phase3_final/coloc/SCIATICA.coloc_variants.tsv.bgz
gzip -cd "$src" | bgzip -@ 4 -c > "$dst.tmp"
mv "$dst.tmp" "$dst"
tabix -f -s 3 -b 4 -e 4 -S 1 "$dst"
md5sum "$dst" > "$dst.md5"
sha256sum "$dst" > "$dst.sha256"
