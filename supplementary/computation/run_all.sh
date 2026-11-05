#!/usr/bin/env bash
# Reproduce all certificate and simulation outputs (except the ECDC analysis, which needs the data file).
set -euo pipefail
cd "$(dirname "$0")"
python scripts/certify_trtu.py
cd scripts
for s in s01_table1_fig1 s02_tables_s1_s2 s03_table_s3_fig2 s04_prop321_plateau s05_resident_first \
         s06_cscp_table2 s07_perturbation s08_table_s4 s09_exact_examples s10_check_thm333_bound; do
  echo "== $s"; python "$s.py"
done
