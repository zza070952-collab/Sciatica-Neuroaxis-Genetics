#!/usr/bin/env bash
set -euo pipefail

project=.
dest="$project/work/software/phase3_twas/FUSION_WEIGHTS_EUR"
mkdir -p "$dest/archives" "$dest/extracted"

manifest="$dest/archive_manifest.tsv"
printf 'tissue\turl\texpected_bytes\tobserved_bytes\tmd5\tsha256\tarchive_test\textracted_pos\n' > "$manifest"

while IFS=$'\t' read -r tissue expected; do
  url="https://s3.us-west-1.amazonaws.com/gtex.v8.fusion/EUR/GTExv8.EUR.${tissue}.tar.gz"
  archive="$dest/archives/GTExv8.EUR.${tissue}.tar.gz"
  curl --fail --location --retry 8 --retry-delay 5 --continue-at - --output "$archive" "$url"
  observed="$(stat -c '%s' "$archive")"
  if [[ "$observed" != "$expected" ]]; then
    printf 'ERROR: %s expected %s bytes, observed %s\n' "$tissue" "$expected" "$observed" >&2
    exit 1
  fi
  tar -tzf "$archive" >/dev/null
  tissue_dest="$dest/extracted/$tissue"
  mkdir -p "$tissue_dest"
  tar -xzf "$archive" -C "$tissue_dest"
  pos_count="$(find "$tissue_dest" -type f -name '*.pos' | wc -l)"
  printf '%s\t%s\t%s\t%s\t%s\t%s\tPASS\t%s\n' \
    "$tissue" "$url" "$expected" "$observed" \
    "$(md5sum "$archive" | awk '{print $1}')" \
    "$(sha256sum "$archive" | awk '{print $1}')" "$pos_count" >> "$manifest"
done <<'EOF'
Nerve_Tibial	251876343
Brain_Spinal_cord_cervical_c-1	208552935
Whole_Blood	195315536
Muscle_Skeletal	207150387
Cells_Cultured_fibroblasts	211305472
EOF

md5sum "$manifest" > "$manifest.md5"
sha256sum "$manifest" > "$manifest.sha256"
