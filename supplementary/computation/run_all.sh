#!/usr/bin/env bash
# Reproduce all certificate and simulation outputs (except the ECDC analysis, which needs the data file).
set -euo pipefail
cd "$(dirname "$0")"
python scripts/certify_trtu.py
python scripts/make_certificate_summary.py
cd scripts
for s in s01_table1_fig1 s02_tables_s1_s2 s03_table_s3_fig2 s04_prop321_plateau s05_resident_first \
         s06_cscp_table2 s07_perturbation s10_check_thm333_bound s11_reimportation \
         s12_return_transmission_finite_start s13_rates_and_shares s14_exact_lemma_A73 \
         s15_revision_checks s16_manuscript_figures s17_finite_variant s18_schematic; do
  echo "== $s"; python "$s.py"
done
