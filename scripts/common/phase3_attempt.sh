#!/usr/bin/env bash
set -u

if [[ $# -lt 8 ]]; then
  echo "usage: $0 MODULE ATTEMPT_ID PRIMARY_OR_FALLBACK DATASET MODEL PARAMETERS OUTPUT_PATH COMMAND..." >&2
  exit 2
fi

module=$1
attempt_id=$2
primary_or_fallback=$3
dataset=$4
model=$5
parameters=$6
output_path=$7
shift 7

project=${PHASE3_PROJECT:-.}
log_tsv="$project/results/phase3_final/00_preflight/phase3_attempt_log.tsv"
log_dir="$project/logs/phase3_final/$module"
mkdir -p "$log_dir" "$output_path" "$(dirname "$log_tsv")"

if [[ ! -e "$log_tsv" ]]; then
  printf 'module\tattempt_id\tprimary_or_fallback\tdataset\tmodel\tparameters\tstart_time\tend_time\tstatus\tfailure_reason\toutput_path\truntime\tpeak_memory\tscientific_usability\n' > "$log_tsv"
fi

start_time=$(date --iso-8601=seconds)
start_epoch=$(date +%s)
stdout="$log_dir/${attempt_id}.stdout.log"
stderr="$log_dir/${attempt_id}.stderr.log"
time_log="$log_dir/${attempt_id}.time.log"
cmd_file="$log_dir/${attempt_id}.command.txt"
printf '%q ' "$@" > "$cmd_file"
printf '\n' >> "$cmd_file"

/usr/bin/time -v -o "$time_log" "$@" >"$stdout" 2>"$stderr"
exit_code=$?
end_epoch=$(date +%s)
end_time=$(date --iso-8601=seconds)
runtime=$((end_epoch-start_epoch))
peak_memory=$(awk -F: '/Maximum resident set size/ {gsub(/^[[:space:]]+/,"",$2); print $2 " kB"}' "$time_log")

if [[ $exit_code -eq 0 ]]; then
  status=SUCCESS
  failure_reason=""
  scientific_usability=PENDING_REVIEW
else
  status=FAILED
  failure_reason=$(tail -n 12 "$stderr" | tr '\t\r\n' ' ' | sed 's/  */ /g; s/^ //; s/ $//')
  scientific_usability=NOT_USABLE
fi

clean_field() { printf '%s' "$1" | tr '\t\r\n' '   '; }
row_file=$(mktemp "$log_dir/${attempt_id}.row.XXXXXX")
{
  clean_field "$module"; printf '\t'
  clean_field "$attempt_id"; printf '\t'
  clean_field "$primary_or_fallback"; printf '\t'
  clean_field "$dataset"; printf '\t'
  clean_field "$model"; printf '\t'
  clean_field "$parameters"; printf '\t'
  clean_field "$start_time"; printf '\t'
  clean_field "$end_time"; printf '\t'
  clean_field "$status"; printf '\t'
  clean_field "$failure_reason"; printf '\t'
  clean_field "$output_path"; printf '\t'
  clean_field "${runtime}s"; printf '\t'
  clean_field "${peak_memory:-NA}"; printf '\t'
  clean_field "$scientific_usability"; printf '\n'
} > "$row_file"
(
  flock -x 9
  cat "$row_file" >> "$log_tsv"
) 9>"$log_tsv.lock"
rm -f "$row_file"
exit "$exit_code"
