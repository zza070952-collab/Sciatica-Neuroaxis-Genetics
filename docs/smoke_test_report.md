# Release-preparation smoke test

Date: 2026-09-05.

Status: PASS (bounded cross-language contract and syntax smoke test).

The curated checkout was transferred to an isolated staging directory on the migrated Linux host. `bash tests/smoke_test.sh` exited 0. Observed output:

```text
{"status": "PASS", "python_syntax_files": 28, "configuration": "PASS", "stage_script_resolution": "PASS", "production_ACAT_identity": "PASS", "production_BH_known_order": "PASS", "full_data_workflow_rerun": false}
R_STARTUP_AND_SYNTAX_PASS
SMOKE_TEST_PASS
```

This verifies Python startup, YAML configuration, all Python production syntax, stage-to-script references, a small execution of the exact production ACAT/BH functions, output-file creation, R startup and R parsing, and Bash syntax. R parsing used the test host's R 4.3.1; the production colocalization record specifies R 4.4.3. This is not an all-environment installation test or a full-data scientific rerun.

The same Python smoke test also passed locally. No research data were simulated and no third-party data were downloaded by these tests.

## Public-disclosure rerun: 2026-09-06

The same retained scientific scripts and test suite passed again on the Linux host: Python syntax/configuration/stage references, production ACAT and BH function checks, R startup and parsing, and Bash syntax all passed. The first attempt used the unconfigured system Python and exited with `ModuleNotFoundError: No module named 'numpy'`. Repeating the test in the installed scientific Python runtime resolved that environment-selection error and returned `SMOKE_TEST_PASS`. Install/activate the documented environment before testing; the operating-system Python is not sufficient. Neither attempt executed a full-data scientific analysis.
