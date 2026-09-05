#!/usr/bin/env python3
"""Stream-audit FinnGen summary statistics and create the common-variant input.

No INFO threshold is applied because the public endpoint file has no INFO field.
Per-variant INFO cannot be reconstructed from this release.
"""

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
from collections import Counter
from datetime import datetime, timezone


def checksum(path, algorithm):
    h = hashlib.new(algorithm)
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--trait", required=True)
    p.add_argument("--output-prefix", required=True)
    p.add_argument("--maf-min", type=float, default=0.01)
    p.add_argument("--palindrome-mid-low", type=float, default=0.40)
    p.add_argument("--palindrome-mid-high", type=float, default=0.60)
    args = p.parse_args()

    os.makedirs(os.path.dirname(args.output_prefix), exist_ok=True)
    out_gz = args.output_prefix + ".common_qc.tsv.gz"
    summary_path = args.output_prefix + ".audit.json"
    counts_path = args.output_prefix + ".qc_counts.tsv"
    start = datetime.now(timezone.utc).astimezone().isoformat()

    counts = Counter()
    chromosomes = Counter()
    header = None
    last_variant = None
    last_position = None
    z_p_n = 0
    z_p_bad = 0
    max_abs_log10_delta = 0.0
    top = []

    with gzip.open(args.input, "rt", newline="") as src, gzip.open(out_gz, "wt", newline="") as dst:
        reader = csv.DictReader(src, delimiter="\t")
        header = reader.fieldnames or []
        required = {"#chrom", "pos", "ref", "alt", "rsids", "pval", "beta", "sebeta", "af_alt"}
        missing = sorted(required.difference(header))
        if missing:
            raise RuntimeError("missing required columns: " + ",".join(missing))
        writer = csv.writer(dst, delimiter="\t", lineterminator="\n")
        writer.writerow(["rsid", "chr", "position", "ref", "effect_allele", "z", "beta", "se", "p", "EAF"])

        for row in reader:
            counts["input_rows"] += 1
            chrom = row["#chrom"]
            chromosomes[chrom] += 1
            ref = row["ref"].upper()
            alt = row["alt"].upper()
            rsid = row["rsids"]
            pos = row["pos"]
            key = (chrom, pos, ref, alt, rsid)
            position_key = (chrom, pos)
            is_exact_duplicate = key == last_variant
            if is_exact_duplicate:
                counts["adjacent_exact_duplicates"] += 1
            if position_key == last_position and key != last_variant:
                counts["adjacent_same_position_other_record"] += 1
            last_variant = key
            last_position = position_key

            try:
                beta = float(row["beta"])
                se = float(row["sebeta"])
                pval = float(row["pval"])
                eaf = float(row["af_alt"])
            except ValueError:
                counts["non_numeric_core"] += 1
                continue
            if not (math.isfinite(beta) and math.isfinite(se) and se > 0 and math.isfinite(pval) and 0 <= pval <= 1 and math.isfinite(eaf) and 0 <= eaf <= 1):
                counts["invalid_core"] += 1
                continue
            z = beta / se
            counts["valid_core"] += 1

            if 1e-250 < pval < 1 and abs(z) < 37:
                p_from_z = math.erfc(abs(z) / math.sqrt(2.0))
                if p_from_z > 0:
                    delta = abs(math.log10(pval) - math.log10(p_from_z))
                    max_abs_log10_delta = max(max_abs_log10_delta, delta)
                    z_p_n += 1
                    if delta > 0.02:
                        z_p_bad += 1

            if len(top) < 100 or pval < top[-1][0]:
                top.append((pval, chrom, int(pos), rsid, beta, se, eaf))
                top.sort(key=lambda x: x[0])
                del top[100:]

            if not rsid.startswith("rs") or any(x in rsid for x in (";", ",", " ")):
                counts["excluded_missing_or_multiple_rsid"] += 1
                continue
            if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT" or ref == alt:
                counts["excluded_non_snv"] += 1
                continue
            maf = min(eaf, 1.0 - eaf)
            if maf < args.maf_min:
                counts["excluded_maf"] += 1
                continue
            if {ref, alt} in ({"A", "T"}, {"C", "G"}) and args.palindrome_mid_low <= eaf <= args.palindrome_mid_high:
                counts["excluded_midfreq_palindrome"] += 1
                continue
            if is_exact_duplicate:
                # Exact duplicates are expected to be adjacent in sorted FinnGen files.
                # The first copy was already written; skip subsequent copies.
                continue
            writer.writerow([rsid, chrom, pos, ref, alt, f"{z:.12g}", row["beta"], row["sebeta"], row["pval"], row["af_alt"]])
            counts["common_qc_rows"] += 1

    end = datetime.now(timezone.utc).astimezone().isoformat()
    summary = {
        "trait": args.trait,
        "input": os.path.abspath(args.input),
        "output": os.path.abspath(out_gz),
        "start_time": start,
        "end_time": end,
        "header": header,
        "info_column_present": any(x.upper() == "INFO" for x in header),
        "info_filter_applied": False,
        "info_gate": "CONDITIONAL_UPSTREAM_QC_ONLY",
        "info_gate_reason": "Public FinnGen endpoint summary statistics omit per-variant INFO; no INFO>=0.8 claim is made.",
        "build": "GRCh38",
        "effect_allele": "ALT",
        "qc_parameters": {
            "maf_min": args.maf_min,
            "snv_only": True,
            "single_rsid_only": True,
            "mid_frequency_palindrome_exclusion": [args.palindrome_mid_low, args.palindrome_mid_high],
        },
        "counts": dict(counts),
        "chromosome_rows": dict(chromosomes),
        "z_p_check": {
            "n_checked": z_p_n,
            "n_abs_log10_delta_gt_0.02": z_p_bad,
            "max_abs_log10_delta": max_abs_log10_delta,
        },
        "top_100_by_p": [
            {"p": x[0], "chr": x[1], "pos": x[2], "rsid": x[3], "beta": x[4], "se": x[5], "eaf": x[6]}
            for x in top
        ],
    }
    with open(summary_path, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    with open(counts_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        for k in sorted(counts):
            writer.writerow([k, counts[k]])
        writer.writerow(["z_p_n_checked", z_p_n])
        writer.writerow(["z_p_bad_delta_gt_0.02", z_p_bad])
        writer.writerow(["z_p_max_abs_log10_delta", max_abs_log10_delta])

    for path in (out_gz, summary_path, counts_path):
        with open(path + ".md5", "w", encoding="ascii") as handle:
            handle.write(f"{checksum(path, 'md5')}  {os.path.basename(path)}\n")
        with open(path + ".sha256", "w", encoding="ascii") as handle:
            handle.write(f"{checksum(path, 'sha256')}  {os.path.basename(path)}\n")


if __name__ == "__main__":
    main()
