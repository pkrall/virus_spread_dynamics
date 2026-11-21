# Simulation code, numerical outputs and interval-arithmetic certificates

This directory contains the code and outputs for

> P. Krall and C. Krall, *Strain elimination and transient dynamics in a clustered population seeded from turnover-renewed transmission units* (submitted to the Journal of Mathematical Biology), and Online Resource 1.

It contains three things:

* **Computer-assisted proofs** (`scripts/certify_trtu.py`). They certify the examples of Lemma A.1.1 (proof in Appendix B.1.1): the fixed points, the convergence and the tabulated extrema for all five scenarios of Table 2 and for all six ramps of Tables S1 and S2. The computation uses interval arithmetic with outward rounding, so every number is an enclosure rather than an estimate.
* **Floating-point computations** (`scripts/s01_*.py` to `scripts/s18_*.py`). They produce the numerical illustrations, Tables 2 and 3, Table S3, all numbers stated in the remarks of Appendix A, and Figs. 1–7 and S1–S3 (Fig. 1 is the model schematic).
* **The ECDC variant analysis** (`scripts/alpha_delta_ordering.py`). It produces Table S6 of Online Resource 1 and its sensitivity analysis.

All outputs under `outputs/` were produced by the code in this directory.

## Requirements and commands

* Python ≥ 3.10. Tested with Python 3.13.18, NumPy 2.5.3, SciPy 1.18.1, Matplotlib 3.11.2 and pandas 3.0.5.
* The certificate script uses only the standard library (`decimal`, `fractions`).
* The other scripts need NumPy, SciPy, Matplotlib and pandas (`requirements.txt`).

```
pip install -r requirements.txt
bash run_all.sh                                   # all outputs except the ECDC analysis; about one minute
python scripts/alpha_delta_ordering.py path/to/strain_data.csv outputs/ecdc
```

`run_all.sh` first runs the certificates (`certify_trtu.py`, then `make_certificate_summary.py`) and then the scripts `s01` to `s18` in order. Each script can also be run on its own from `scripts/`.

## Model conventions

* **Transmission unit (TRTU).** One update step per time step: p_S(t+1) = (1 − m)(p_S(t) − p_N(t)) + m, p_I(t) = Σ_{j = D_E+1}^{D_E+D_I} (1 − m)^j p_N(t − j), p_N(t) = p_S(t)[1 − (1 − τ p_I(t))^K]. With two strains, the exposure probability is allocated in proportion to τ_V p_{I,V}(t).
* **Seeding.** A naive unit is seeded with 0.01 m at step 0. In replacement experiments, Y is seeded with 0.01 m into a unit whose strain X is at its fixed point (constant fixed-point history before step 0).
* **Delayed state.** X_t = (p_S(t), p_N(t−1), …, p_N(t−L)) with L = D_E + D_I, and X_{t+1} = F(X_t).
* **Timing conventions (Tables 3 and S3).** *Counted* steps include both the first step at which the Y share reaches 1 % and the first step at which it reaches 95 %. *Elapsed* steps are differences of step indices.
* **Seeded population (CSCP).** Expected response kernels (Section 2.3) are computed exactly by propagating the full distribution of the within-cluster chain binomial; cross-cluster seeding is added by a renewal convolution. The exposure intensity is the convolution of the seeding with the response kernel (Corollary A.2.3). Stationary prehistories are represented by a constant seeding prefix longer than the support of the kernel. Exposure intensities are per individual in the transmission unit (coupling coefficient α = 0.1, population of the transmission unit normalized to one).

## Article item → script → output

| Article / Online Resource 1 | Script | Output (`outputs/…`) |
|---|---|---|
| Lemma A.1.1, Appendix B.1.1: certified fixed points, contraction boxes and extrema for Table 2 | `certify_trtu.py` | `certificates/fixed_points.csv`, `contraction.csv`, `table1_certificate.csv`, `tube_scenario{1..5}.csv`, `enclosures_table1_scenario{1..5}.csv`, `summary.txt` |
| Tables S1, S2 (ramps of K), certified | `certify_trtu.py` | `certificates/tables_S1_S2_certificate.csv`, `enclosures_tableS{1,2}_ramp{5,10,20}.csv` |
| Certificate summary (box radius, iterate power, derivative bound, entry step) | `make_certificate_summary.py` | `certificates/certificate_summary.csv` |
| Table 2 (floating point, with step indices) | `s01_table1_fig1.py` | `tables/table1.csv` |
| Tables S1, S2 (floating point), Figs. S1, S2 | `s02_tables_s1_s2.py` | `tables/tables_S1_S2.csv`, `figures/figS1.png`, `figS2.png` |
| Table S3; replacement times and total dips quoted with Fig. 3 and in the remark after Corollary A.4.1 | `s03_table_s3_fig2.py` | `tables/table_S3.csv` |
| Lemma A.6.1, numerical illustration (plateauing), Fig. 6 | `s04_prop321_plateau.py` | `tables/prop321_plateau_example.csv`, `figures/fig6.png` |
| Plateau during strain replacement (after Lemma A.6.1), Fig. S3 | `s05_resident_first.py` | `tables/prop321_resident_first_example.csv`, `figures/figS3.png` |
| Table 3 (in-cluster probabilities 0.06094531916, 0.08131295535) | `s06_cscp_table2.py` | `tables/table2.csv`, `table2_profiles.csv` |
| Corollary A.7.2, numerical illustration (fading of perturbations) | `s07_perturbation.py` | `tables/cor336_perturbation.csv` |
| Lemma A.2.2(a): numerical check of the bound 2A/B (not part of a proof) | `s10_check_thm333_bound.py` | `tables/thm333a_numerical_check.csv` |
| Remark (re-importation of the incumbent) after Lemma A.1.2 | `s11_reimportation.py` | `tables/remark_reimportation.csv` |
| Return-transmission experiment (remark after Lemma A.1.2); finite-start accounting (remark after Lemma A.3.1) | `s12_return_transmission_finite_start.py` | `tables/remark_return_transmission.csv`, `remark_finite_start.csv` |
| Remark on the rate of elimination (0.906, 0.883); Corollary A.4.1(f) (decay factors 0.981, 0.963, 0.969); Lemma A.2.2(a) (0.175 and 0.145); shares after Lemma A.3.1 (70 %, 78 %, 82 %, 91 %); replacement-window accounting (67–68 %, 15 exposures, 2.4 and 3.5 clusters per individual in the transmission units) | `s13_rates_and_shares.py` | `tables/rates_and_shares.csv`, `replacement_window_accounting.csv` |
| Proof of Lemma A.7.3 (Appendix B.2.8): exact check of the three-system construction | `s14_exact_lemma_A73.py` | `tables/lemma_A73_exact.csv` |
| Checks for the revision (no numbers printed in the article): Section 3.2(ii) conservation, remark on partially immune arrivals after Theorem 5, Theorem 4(e) | `s15_revision_checks.py` | `tables/revision_checks.csv` |
| Table 4 and the remark (scale of the approximation) after Theorem 5: finite-population variant of Section 3.2 for the finite-start example (X at step 0, Y at step 400; r = 0, 2·10⁻⁵, 2·10⁻⁴, 10⁻³ and 1/(100n), that is, 100 clusters per individual in the transmission units) | `s17_finite_variant.py` | `tables/table4_theorem5_scale.csv` |
| Figs. 2, 3, 4, 5 and 7 | `s16_manuscript_figures.py` | `figures/fig2.png`, `fig3.png`, `fig4.png`, `fig5.png`, `fig7.png` |
| Fig. 1 (model schematic; drawing only) | `s18_schematic.py` | `figures/fig1.png` |
| Table S6 and sensitivity analysis (Online Resource 1, Section S2.2) | `alpha_delta_ordering.py` | `ecdc/table_s6.csv`, `ecdc/sensitivity.csv` |

Some file names keep an earlier numbering (`s01_table1_fig1.py` and `s03_table_s3_fig2.py` relate to Table 2 and Fig. 2, and to Table S3 and Fig. 3; `table1.csv` is Table 2 of the article, `table2.csv` is Table 3, `prop321_*` refers to Lemma A.6.1, `cor336_*` to Corollary A.7.2, `thm333a_*` to Lemma A.2.2(a)). Files named `figures/work_*.png` are working plots and do not appear in the article.

## The computer-assisted proofs (`certify_trtu.py`)

### Arithmetic

* **Representation and rounding.** `trtu_model/interval.py` implements closed intervals with `decimal.Decimal` endpoints at 80 significant digits. Lower endpoints are computed with `ROUND_FLOOR` and upper endpoints with `ROUND_CEILING`, so every elementary operation (+, −, ×, ÷, integer powers) returns an interval that contains the exact real result.
* **Parameters.** Parameters are entered as decimal strings, for example τ = `'0.12'` and m = `'0.035'`, and are therefore exact rather than binary approximations.
* **Printing.** Printed bounds are rounded outward again, so a printed interval always contains the computed one; the derivative bounds in `certificate_summary.csv` are rounded upward.

### Steps of each certificate

1. **Fixed point.** A floating-point estimate is refined to a bracket [a, b] of width about 10⁻⁴⁰. Interval evaluation proves H(a) > 0 > H(b) for the fixed-point equation H, so the intermediate value theorem gives a fixed point p_N* ∈ [a, b] (`fixed_points.csv`). The bracket printed in Appendix B.1.1 is checked separately (`summary.txt`).
2. **Contraction box.** Let c be the delayed state at the fixed point and B the closed sup-norm box of radius r around c. The script encloses the images B_k ⊇ F^k(B), k = 0…P, and multiplies interval enclosures of the Jacobians DF over B_{P−1}, …, B_0. This proves ‖DF^P(x)‖_∞ ≤ L < 1 for all x ∈ B. It also checks L·(r + w) + w ≤ r, where w bounds the distance between c and the exact fixed point; together these give F^P(B) ⊆ B. Every trajectory that enters B therefore stays in the tube ⋃_k B_k and converges to the fixed point (`contraction.csv`, `tube_scenario*.csv`).
3. **Trajectory enclosure.** The trajectory from the seeding (Table 2), or from the old fixed point under the ramp of K (Tables S1, S2), is enclosed step by step until it provably lies in B (`enclosures_*.csv`).
4. **Extrema.** Each tabulated extremum is certified by strict separation: its enclosure lies strictly above (or below) the enclosures of all competing steps and the range of p_N over the tube, which covers all later steps. Rounding the enclosure to four decimals must give the tabulated value; for Tables S1 and S2 the step index must also match (`table1_certificate.csv`, `tables_S1_S2_certificate.csv`; the column `certified_and_matches` is `True` in every row).

### Certificate summary (`certificates/certificate_summary.csv`)

| Certified trajectory | box of scenario | box radius r | P | bound on ‖DF^P‖_∞ over B | box entered at step |
|---|---|---|---|---|---|
| Table 2, scenario 1 | 1 | 1·10⁻⁵ | 78 | 0.9072 | 588 |
| Table 2, scenario 2 | 2 | 5·10⁻⁵ | 27 | 0.8869 | 229 |
| Table 2, scenario 3 | 3 | 5·10⁻⁵ | 24 | 0.8462 | 203 |
| Table 2, scenario 4 | 4 | 5·10⁻⁵ | 15 | 0.8956 | 123 |
| Table 2, scenario 5 | 5 | 5·10⁻⁵ | 21 | 0.8917 | 154 |
| Table S1, ramps of 5, 10, 20 steps | 4 | 5·10⁻⁵ | 15 | 0.8956 | 92, 95, 99 |
| Table S2, ramps of 5, 10, 20 steps | 3 | 5·10⁻⁵ | 24 | 0.8462 | 161, 163, 148 |

The full certificate run takes a few seconds.

## Floating-point computations

* **Precision.** Double precision. Where an item is also certified (Table 2, Tables S1–S2), the floating-point and certified values agree to all printed digits.
* **Experiments outside the model.** `s11` (re-importation) uses the resident seeding rule of Section 2.1 at every step; `s12` (return transmission) makes exposed arrivals proportional to the exposure intensity of the seeded population, σ_V(t) = φ ℐ_V(t−1) with φ = 0.001, 0.01, 0.1. Neither is part of the model for which the theorems are proved. `s17` solves the finite-population variant of Section 3.2, to which Theorem 5 refers.
* **Perturbation illustration.** The values printed in the article correspond to constant seeding κ = 0.025 (baseline exposure intensity 0.25). The expected dynamics are linear in κ, so ratios and step counts are the same at every seeding level.

## ECDC data (Table S6)

* **Input not included.** The ECDC file "Data on SARS-CoV-2 variants in the EU/EEA" is not redistributed here.
* **File used.** The analysis used a download with SHA-256 `f20a1f1cdd34b2ac171ae42ed1bb6ccd54f8e31bfe283eeaf306250d9f4acf4d`, semicolon-separated as saved from a spreadsheet program. The script also accepts the comma-separated original.
* **Definitions.** Lineage definitions, thresholds and the endpoint rule are documented at the top of the script.

## Licence and citation

[to be added by the authors]
