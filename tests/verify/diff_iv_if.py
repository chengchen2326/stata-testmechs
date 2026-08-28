#!/usr/bin/env python3
"""Compare IV IF from R and Stata dumps.

Requires prior runs of:
  Rscript tests/verify/dump_iv_if_R_all.R
  stata -e "do tests/verify/dump_iv_if_stata_all.do"

Outputs pass/fail per instrument.
"""
import csv
import sys

TOL = 1e-6  # tolerance: matches must be within 1e-6 to PASS

exit_code = 0
for instr in ["treat2", "treat3", "treat4"]:
    try:
        with open(f"/tmp/r_iv_if_{instr}.csv") as f:
            r_if = [float(row["if_treat"]) for row in csv.DictReader(f)]
        with open(f"/tmp/stata_iv_if_{instr}.csv") as f:
            s_if = [float(row["if_treat1"]) for row in csv.DictReader(f)]
    except FileNotFoundError as e:
        print(f"{instr}: MISSING FILE: {e.filename}")
        exit_code = 1
        continue

    if len(r_if) != len(s_if):
        print(f"{instr}: LEN MISMATCH r={len(r_if)} s={len(s_if)}")
        exit_code = 1
        continue

    diffs = [abs(a - b) for a, b in zip(r_if, s_if)]
    max_d = max(diffs)
    mean_d = sum(diffs) / len(diffs)
    if max_d < TOL:
        status = "PASS"
    elif max_d < 1e-4:
        status = "CLOSE"
        exit_code = 1
    else:
        status = "FAIL"
        exit_code = 1
    print(f"{instr}: max|diff|={max_d:.4e}  mean|diff|={mean_d:.4e}  [{status}]")

sys.exit(exit_code)
