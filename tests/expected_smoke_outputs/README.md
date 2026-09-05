The smoke run writes `results/smoke/smoke_output.json` with status PASS and checks
configuration, script references, Python parsing and the production ACAT/BH functions.
`smoke_test.sh` also starts R, parses every R script and checks shell syntax.
Expected analytical unit checks: ACAT(0.5, 0.5) = 0.5;
BH(0.01, 0.04, 0.03) = (0.03, 0.04, 0.04).
These are software tests, not study observations. The smoke test does not claim
to rerun genome-wide analyses or validate all installed scientific dependencies.
