# TRTU / TRTA / CSCP model: simulation code and certificate outputs

This repository contains code and output tables for the following article and its online resource:

> P. Krall and C. Krall, *Strain replacement independent of susceptibility depletion in subcritical clustered populations seeded by turnover-renewed transmission units* (submitted to the Journal of Mathematical Biology), and Online Resource 1.

It contains three things:

* **Computer-assisted proofs** (`scripts/certify_trtu.py`). These cover Theorem 3.1.3 and Appendix A.1.3, for all five scenarios of Table 1 and all six ramps of Tables S1 and S2. The computation uses rigorous interval arithmetic with outward rounding, so every number is an enclosure rather than an estimate.
* **Floating-point simulations** (`scripts/s01_*.py` to `scripts/s10_*.py`). These produce the numerical illustrations, Tables 2 and S3–S4, and Figs. 1–6 and S2–S4.
* **The ECDC variant analysis** (`scripts/alpha_delta_ordering.py`). This produces Table S7 and its sensitivity analysis.

All numbers below were produced by the code in this repository; the outputs are included under `outputs/`.

## Requirements

* Python ≥ 3.10.
* The certificate script uses only the standard library: `decimal` and `fractions`.
* The simulations and the ECDC analysis need NumPy, SciPy, Matplotlib and pandas (see `requirements.txt`).
* Tested versions: Python 3.13.18, NumPy 2.5.3, SciPy 1.18.1, Matplotlib 3.11.2, pandas 3.0.5.

```
pip install -r requirements.txt
bash run_all.sh                     # everything except the ECDC analysis (Table S4's Monte Carlo takes a few minutes)
python scripts/alpha_delta_ordering.py path/to/strain_data.csv outputs/ecdc
```

## Model conventions used throughout

* **TRTU update.** The time step is one update step. The recursion is p_S(t+1) = (1 − m)(p_S(t) − p_N(t)) + m, with:
  * p_I(t) = Σ_{j = D_E+1}^{D_E+D_I} (1 − m)^j p_N(t − j);
  * p_N(t) = p_S(t)[1 − (1 − τ p_I(t))^K].
* **Seeding and steps.** A naive unit is seeded with 0.01 m at step 0, so p_N(0) = 0.01 m and p_S(0) = 1. Steps are counted from the seeding.
* **Delayed state.** The state is X_t = (p_S(t), p_N(t−1), …, p_N(t−L)) with L = D_E + D_I, and X_{t+1} = F(X_t).
* **Replacement times (Tables 2 and S3).** Y is seeded with 0.01 m at step 50 into a unit whose strain X is at its fixed point.
  * *Counted* steps include both the first step at which the Y share reaches 1 % and the first step at which it reaches 95 %.
  * *Elapsed* steps are differences of step indices.
* **CSCP profiles.**
  * Expected cascade-age profiles are computed exactly by propagating the full distribution of the within-cluster chain binomial.
  * Cross-cluster seeding is added by a renewal convolution.
  * The CSCP exposure intensity is the convolution of the seeding with this profile (Corollary 3.3.4).
  * Stationary prehistories (shift extension) are represented by a constant seeding prefix longer than the support of the profile.

## Paper item → script → output

| Article / Online Resource 1 | Script | Output (`outputs/…`) |
|---|---|---|
| Theorem 3.1.3, Appendix A.1.3: fixed points, contraction boxes, certified extrema for Table 1 | `certify_trtu.py` | `certificates/fixed_points.csv`, `contraction.csv`, `table1_certificate.csv`, `tube_scenario{1..5}.csv`, `enclosures_table1_scenario{1..5}.csv`, `summary.txt` |
| Tables S1, S2 (linear ramps of K), certified | `certify_trtu.py` | `certificates/tables_S1_S2_certificate.csv`, `enclosures_tableS{1,2}_ramp{5,10,20}.csv` |
| Table 1, Fig. 1 (floating point, with step indices) | `s01_table1_fig1.py` | `tables/table1.csv`, `figures/fig1_*.png` |
| Tables S1, S2, Figs. S2, S3 (floating point) | `s02_tables_s1_s2.py` | `tables/tables_S1_S2.csv`, `figures/figS2_*.png`, `figS3_*.png` |
| Table S3, Fig. 2 | `s03_table_s3_fig2.py` | `tables/table_S3.csv`, `figures/fig2_*.png` |
| Section 3.2.1 plateau illustration, Fig. 3 | `s04_prop321_plateau.py` | `tables/prop321_plateau_example.csv`, `figures/fig3_*.png` |
| Resident-first plateau during replacement, Fig. S4 | `s05_resident_first.py` | `tables/prop321_resident_first_example.csv`, `figures/figS4_*.png` |
| Table 2, Figs. 4, 5 (in-cluster probabilities 0.06094531916, 0.08131295535) | `s06_cscp_table2.py` | `tables/table2.csv`, `table2_profiles.csv`, `figures/fig4_*.png`, `fig5_*.png` |
| Corollary 3.3.6 illustration, Fig. 6 | `s07_perturbation.py` | `tables/cor336_perturbation.csv`, `figures/fig6_*.png` |
| Table S4 (Corollary 3.3.8) | `s08_table_s4.py` | `tables/table_S4.csv`, `table_S4_cascade_sizes.csv` |
| Proposition 3.3.7(a), exact example | `s09_exact_examples.py` | `tables/prop337a_exact.csv` |
| Theorem 3.3.3(a): numerical check of the bound (not part of the proof) | `s10_check_thm333_bound.py` | `tables/thm333a_numerical_check.csv` |
| Table S7, sensitivity analysis (Section S4) | `alpha_delta_ordering.py` | `ecdc/table_s7.csv`, `ecdc/sensitivity.csv` |

## The computer-assisted proofs (`certify_trtu.py`)

### Arithmetic

* **Representation and rounding.** `trtu_model/interval.py` implements closed intervals with `decimal.Decimal` endpoints at 80 significant digits.
  * Lower endpoints are computed with `ROUND_FLOOR` and upper endpoints with `ROUND_CEILING`.
  * Every elementary operation (+, −, ×, ÷, integer powers) therefore returns an interval that contains the exact real result.
* **Parameters.** The parameters are entered as decimal strings, e.g. τ = `'0.12'` and m = `'0.035'`. They are therefore exact and not binary approximations.
* **Printing.** Printed bounds are rounded outward again, so a printed interval always contains the computed one.

### Steps of each certificate

1. **Fixed point.**
   * A floating-point estimate is refined to a bracket [a, b] of width about 10⁻⁴⁰.
   * Interval evaluation proves H(a) > 0 > H(b) for the fixed-point equation H, so the intermediate value theorem gives a fixed point p_N* ∈ [a, b].
   * The output is `fixed_points.csv`. The bracket printed in Appendix A.1.3 is checked separately; see `summary.txt`.
2. **Contraction box.**
   * Let c be the delayed state at the fixed point, and let B be the closed sup-norm box of radius r around c.
   * The script encloses the images B_k ⊇ F^k(B), k = 0…P, and multiplies interval enclosures of the Jacobians DF over B_{P−1}, …, B_0.
   * This proves ‖DF^P(x)‖_∞ ≤ L < 1 for all x ∈ B. The script also checks L·(r + w) + w ≤ r, where w bounds the distance between c and the exact fixed point. Together these give F^P(B) ⊆ B.
   * Every trajectory that enters B therefore stays in the tube ⋃_k B_k and converges to the fixed point.
   * The output is `contraction.csv`, with the full tube in `tube_scenario*.csv`.
3. **Trajectory enclosure.**
   * The trajectory from the seeding (Table 1), or from the old fixed point under the K ramp (Tables S1, S2), is enclosed step by step until the state provably lies in B.
   * Every step is recorded in `enclosures_*.csv`.
4. **Extrema.**
   * Each tabulated extremum is certified by strict separation. Its enclosure lies strictly above (or below) the enclosures of all competing steps, and also strictly above (or below) the range of p_N over the tube, which covers all later steps.
   * Rounding the enclosure to four decimals must give the tabulated value. For Tables S1 and S2 the step index must also match.
   * The outputs are `table1_certificate.csv` and `tables_S1_S2_certificate.csv`. The column `certified_and_matches` is `True` in every row.

### Contraction parameters found

| Scenario | box radius r | P | bound on ‖DF^P‖_∞ over B | box entered (Table 1) |
|---|---|---|---|---|
| 1 | 1·10⁻⁵ | 78 | 0.9072 | step 588 |
| 2 | 5·10⁻⁵ | 27 | 0.8869 | step 229 |
| 3 | 5·10⁻⁵ | 24 | 0.8462 | step 203 |
| 4 | 5·10⁻⁵ | 15 | 0.8956 | step 123 |
| 5 | 5·10⁻⁵ | 21 | 0.8917 | step 154 |

The ramps of Table S1 end in scenario 4's box, and those of Table S2 in scenario 3's box. The full certificate run takes a few seconds.

## Floating-point simulations

* **Precision.** The simulations use double precision. Where an item is also certified (Table 1, Tables S1–S2), the floating-point and certified values agree to all printed digits.
* **Table S4, left group: Monte Carlo.** This group has n = 10, D_I = 20 and perfect isolation. The exact within-cluster state space is too large there, so the cascade-age profile is estimated by Monte Carlo:
  * first-passage percolation with geometric transmission delays;
  * 2,000,000 clusters per strain;
  * NumPy generator `default_rng(1)`.

  The expected final sizes of this group are also given exactly (Reed–Frost) in `table_S4_cascade_sizes.csv`.
* **Table S4, right group: exact.**
* **Corollary 3.3.6 illustration.** The values printed in the article correspond to constant seeding κ = 0.025, which gives a baseline exposure intensity of 0.25. The expected field is linear in κ, so peak/baseline ratios and step counts are the same at every seeding level; `cor336_perturbation.csv` lists both.

## ECDC data (Table S7)

* **Input not included.** The ECDC file "Data on SARS-CoV-2 variants in the EU/EEA" is not redistributed here.
* **File used.** The analysis used a download with SHA-256 `f20a1f1cdd34b2ac171ae42ed1bb6ccd54f8e31bfe283eeaf306250d9f4acf4d`, semicolon-separated as saved from a spreadsheet program. The script also accepts the comma-separated original.
* **Definitions.** Lineage definitions, thresholds and the endpoint rule are documented at the top of the script.

## Licence and citation

[to be added by the authors]
