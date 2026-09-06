# Security and privacy audit

Status: PASS

- Reviewed text files: 92
- Largest file: `scripts/11_figures_computational/main_figures.py` (45454 bytes)
- Files above 50 MiB: 0
- Credential patterns, secret assignments, personal absolute paths, phone/email patterns, binary/data payloads and internal script references were checked.
- Git history was empty at initial remote inspection; only this curated tree is eligible for the initial commit.
- No wet-laboratory data or statistical scripts, third-party matrices/weights/LD, cache or large log files are included.
- README descriptions of excluded data are documentation, not included data.

## Public-disclosure recheck: 2026-09-06

Status: PASS. The tracked tree and 93 unique blobs from the existing two-commit history were scanned again. No credential, secret assignment, personal path, phone, email, IPv4 literal, prohibited binary/data payload, or file above 50 MiB was detected. Public disclosure was explicitly authorized by the author. This is a security review, not a certification of complete scientific reproducibility; the latter limitations remain in the README and reproducibility documentation.
