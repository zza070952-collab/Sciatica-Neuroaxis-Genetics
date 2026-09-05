#!/usr/bin/env python3
import csv
import hashlib
import re
from pathlib import Path

ROOT = Path(".")
OUT = ROOT / "results/phase3_final/02_genetic_architecture"
LOG = OUT / "ldsc_logs"

def grab(pattern, text, flags=0):
    m = re.search(pattern, text, flags)
    return m.groups() if m else None

rows = []
for path in sorted(LOG.glob("*_observed.log")) + sorted(LOG.glob("*_liability_*.log")):
    text = path.read_text(errors="replace")
    trait = path.name.split("_")[0]
    scale = "observed" if "_observed" in path.name else "liability_sensitivity"
    prev = "NA"
    if scale.startswith("liability"):
        prev = {"K003": "0.03", "K005": "0.05", "K010": "0.10"}[path.stem.rsplit("_", 1)[-1]]
    obs = grab(r"Total Observed scale h2: ([^ ]+) \(([^)]+)\)", text)
    lia = grab(r"Total Liability scale h2: ([^ ]+) \(([^)]+)\)", text)
    lam = grab(r"Lambda GC: ([^\n]+)", text)
    mean = grab(r"Mean Chi\^2: ([^\n]+)", text)
    intercept = grab(r"Intercept: ([^ ]+) \(([^)]+)\)", text)
    ratio = grab(r"Ratio: ([^ ]+) \(([^)]+)\)", text)
    if not (obs or lia) or not intercept:
        raise RuntimeError(f"failed to parse {path}")
    h, se = lia if lia else obs
    rows.append([trait, scale, prev, h, se, float(h)/float(se), lam[0], mean[0], *intercept, *(ratio or ("NA", "NA")), str(path.relative_to(ROOT))])

with open(OUT / "ldsc_heritability.tsv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["trait", "scale", "population_prevalence_assumption", "h2", "se", "z", "lambda_gc", "mean_chisq", "intercept", "intercept_se", "ratio", "ratio_se", "evidence_log"])
    w.writerows(rows)

rg_text = (LOG / "SCIATICA_LDH_rg.log").read_text(errors="replace")
rg = grab(r"Genetic Correlation: ([^ ]+) \(([^)]+)\)\s+Z-score: ([^\n]+)\s+P: ([^\n]+)", rg_text)
gcov = grab(r"Total Observed scale gencov: ([^ ]+) \(([^)]+)\)", rg_text)
gcint = grab(r"Genetic Covariance\s+-+.*?Intercept: ([^ ]+) \(([^)]+)\)", rg_text, re.S)
if not rg or not gcov or not gcint:
    raise RuntimeError("failed to parse rg log")
with open(OUT / "ldsc_genetic_correlation.tsv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["trait1", "trait2", "rg", "se", "z", "p", "genetic_covariance", "genetic_covariance_se", "sampling_covariance_intercept", "sampling_covariance_intercept_se", "evidence_log"])
    w.writerow(["SCIATICA", "LDH", *rg, *gcov, *gcint, str((LOG / "SCIATICA_LDH_rg.log").relative_to(ROOT))])

(OUT / "genetic_architecture_report.md").write_text("""# Genetic architecture of Sciatica and LDH

Primary LDSC estimates are reported on the observed scale using European HapMap3 summary statistics and European 1000 Genomes LD scores. Binary-trait liability estimates are sensitivity analyses because no single population prevalence was prespecified: Sciatica uses K=0.05 and 0.10; LDH uses K=0.03 and 0.05. These assumptions change h2 scale but not the genetic-correlation interpretation.

The bivariate run explicitly records the genetic covariance and the cross-trait LDSC intercept. Because the two FinnGen endpoints share release participants and controls, the cross-trait intercept is interpreted as sampling covariance/sample overlap rather than independent replication.

No three-trait GenomicSEM or Q-SNP analysis was run in Phase 3.
""", encoding="utf-8")

def digest(path, name):
    h = hashlib.new(name)
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(8 * 1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()

for p in OUT.rglob("*"):
    if p.is_file() and p.suffix not in {".md5", ".sha256"}:
        Path(str(p)+".md5").write_text(f"{digest(p,'md5')}  {p.name}\n")
        Path(str(p)+".sha256").write_text(f"{digest(p,'sha256')}  {p.name}\n")
