#!/usr/bin/env bash
set -euo pipefail
root=.
cd "$root"
mkdir -p work/phase3_final/coloc
src=work/gwas_qc/SCIATICA_R13.common_qc.tsv.gz
dst=work/phase3_final/coloc/SCIATICA.coloc_eaf.tsv.bgz
gzip -cd "$src" | bgzip -@ 4 -c > "$dst.tmp"
mv "$dst.tmp" "$dst"
tabix -f -s 2 -b 3 -e 3 -S 1 "$dst"
md5sum "$src" "$dst" > work/phase3_final/coloc/SCIATICA.coloc_eaf.source_and_indexed.md5
sha256sum "$src" "$dst" > work/phase3_final/coloc/SCIATICA.coloc_eaf.source_and_indexed.sha256
