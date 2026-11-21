"""Checks for statements added in the revision (no numbers from this script are printed in the article):
(a) Section 3.2(ii): with infected and immune arrivals the compartment fractions of a transmission unit sum to one and
    stay nonnegative (single strain, scenario 1, constant a and z);
(b) remark after Theorem 5: with a constant immune fraction z among arrivals, a unit is endemic iff (1-z) R_df > 1,
    and the endemic level decreases in z (simulation vs fixed-point formula);
(c) Theorem 4(e): cluster-level reproduction number 1 - Gamma_C/Gamma at Gamma = 10, and its value after the
    cross-cluster rate is multiplied by Gamma/(Gamma-1)."""
import numpy as np
from _common import *

rows = []
# (a)
tau, K, m, DE, DI = SCEN[1]; q = 1 - m; L = DE + DI
for a, z in ((0.0, 0.0), (0.05, 0.1), (0.2, 0.3)):
    S = 1.0; R = 0.0; ages = np.zeros(L + 1); worst = 0.0; mn = 1.0
    for t in range(2000):
        I = ages[DE + 1:L + 1].sum()
        N = S * (1 - (1 - tau * I) ** K) + (0.01 * m if t == 5 else 0.0)
        S = (S - N) * q + m * (1 - a - z)
        R = (R + ages[L]) * q + m * z
        ages[1:] = ages[:-1] * q; ages[1] = (N + m * a / q) * q
        tot = S + ages[1:].sum() + R; worst = max(worst, abs(tot - 1)); mn = min(mn, S, R, ages.min())
    rows.append([f'(a) arrivals exposed {a}, immune {z}: max |sum of compartments - 1|', f'{worst:.1e}'])
    rows.append([f'(a) arrivals exposed {a}, immune {z}: minimum compartment', f'{mn:.3g}'])
# (b)
for sc in (1, 3):
    tau, K, m, DE, DI = SCEN[sc]; q = 1 - m; d = DE + 1; L = DE + DI
    AI = sum(q ** j for j in range(d, L + 1)); Rdf = tau * K * AI; zc = 1 - 1 / Rdf
    rows.append([f'(b) scenario {sc}: R_df, threshold 1 - 1/R_df', f'{Rdf:.4f}, {zc:.4f}'])
    for z in (0.0, 0.1, zc - 0.05, zc - 0.01, zc + 0.02):
        off = L + 2; N = np.zeros(40000 + off); S = np.ones(40000 + off) * (1 - z); N[:off] = 1e-3
        for t in range(off, 40000 + off):
            S[t] = (S[t - 1] - N[t - 1]) * q + m * (1 - z)
            N[t] = S[t] * (1 - (1 - tau * sum(q ** j * N[t - j] for j in range(d, L + 1))) ** K)
        f = lambda x: ((1 - z) - q / m * x) * (1 - (1 - tau * AI * x) ** K) - x
        lo, hi = 1e-14, m
        if f(lo) > 0:
            for _ in range(200):
                mid = (lo + hi) / 2; lo, hi = (mid, hi) if f(mid) > 0 else (lo, mid)
            fp = lo
        else: fp = 0.0
        rows.append([f'(b) scenario {sc}, z = {z:.3f}, (1-z)R_df = {(1 - z) * Rdf:.3f}: simulated p_N, fixed point', f'{N[-1]:.3e}, {fp:.3e}'])
# (c)
G = 10.0
for name, GC in (('x1 (mass action)', 1.0), ('x2', 4.0), ('x3', 3.0)):
    Rc = 1 - GC / G; Rc2 = Rc * G / (G - 1)
    rows.append([f'(c) {name}: R_c at Gamma = 10; after rho x Gamma/(Gamma-1)', f'{Rc:.4f}; {Rc2:.4f}'])
write_csv('revision_checks.csv', ['check', 'value'], rows)
