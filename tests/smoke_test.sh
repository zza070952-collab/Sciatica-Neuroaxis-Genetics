#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
python3 tests/test_smoke.py
python3 workflow/run.py --dry-run > results/smoke/workflow_contract.json
Rscript -e 'stopifnot(getRversion() >= "4.3.0"); x <- list.files("scripts", pattern="[.]R$", recursive=TRUE, full.names=TRUE); for (p in x) parse(p); cat("R_STARTUP_AND_SYNTAX_PASS\n")'
while IFS= read -r f; do bash -n "$f"; done < <(find scripts -name '*.sh' -type f)
printf 'SMOKE_TEST_PASS\n'
