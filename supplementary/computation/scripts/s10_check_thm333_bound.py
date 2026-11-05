"""Numerical check (not part of a proof) of the bound in Theorem 3.3.3(a): ratio of the limiting decay rates
against 2A/B over a parameter grid, and depletion-length ratios along a family with in-cluster probability -> 1."""
import numpy as np
from scipy.optimize import brentq
from _common import *
from trtu_model.cscp import total_profile, cluster_profile
np.seterr(over='ignore')   # brentq probes large u; overflow to inf is harmless there
M = lambda u, DE, DI: np.mean([np.exp(u * j) for j in range(DE + 1, DE + DI + 1)])
worst = 0
for G in (2.5, 3, 5, 10, 30, 100, 1000):
    for n in range(2, 12):
        if not n < G: continue
        for DE in (1, 2, 3, 5, 10):
            for DI in (2, 3, 4, 7, 15, 30):
                A = np.log(G / (G - 1)); B = np.log(G / (G - n))
                mu1 = brentq(lambda u: M(u, DE, DI) - np.exp(A), 1e-14, 50)
                muinf = brentq(lambda u: (1 + (n - 1) * np.exp(u * (DE + 1))) / n * M(u, DE, DI) - np.exp(B), 1e-14, 50)
                worst = max(worst, (mu1 / muinf) / (2 * A / B))
G, n, DE, DI = 10.0, 5, 3, 4
def depletion(P, f):
    tail = P['Gamma_tot'] - P['Gam']; return int(np.where(tail / P['Gamma_tot'] <= f)[0][0]) + 1
P1 = total_profile(1, 0.0, (G - 1) / (G * DI), DE, DI, A=4000); rows = [['max over grid of (mu_1/mu_inf)/(2A/B)', round(worst, 6), '']]
for rb in (0.5, 0.9, 0.99, 0.999):
    GC = cluster_profile(n, rb, DE, DI)[0][-1]; rho = (G - GC) / (G * GC * DI)
    Pi = total_profile(n, rb, rho, DE, DI, A=4000)
    rows.append([f'in-cluster probability {rb}: depletion-length ratio for f = 1e-3, 1e-6, 1e-9, 1e-12', [round(depletion(Pi, f) / depletion(P1, f), 4) for f in (1e-3, 1e-6, 1e-9, 1e-12)], round(Pi['Gamma_tot'], 6)])
write_csv('thm333a_numerical_check.csv', ['check', 'value', 'Gamma_total'], rows)
