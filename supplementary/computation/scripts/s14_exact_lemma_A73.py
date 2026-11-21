"""Exact check of the three-system construction in the proof of Lemma A.7.3 (Appendix B.2.8):
D_E = 1, D_I = 2, multiplier f = 2, G3(p) = 1 + 2(q + q^2 - q^3) with q = 1 - (1 - p)^2 for a three-member cluster.
Rational arithmetic for the inequalities; the within-cluster chain is used as a floating-point cross-check of G3."""
from fractions import Fraction as F
from math import sqrt
from _common import *
from trtu_model.cscp import cluster_profile

DE, DI = 1, 2
G3 = lambda p: 1 + 2 * ((1 - (1 - p) ** 2) + (1 - (1 - p) ** 2) ** 2 - (1 - (1 - p) ** 2) ** 3)
rows = []
for p in (F(11, 100), F(3, 25), F(11, 50)):
    rows.append([f'G3({p})', str(G3(p)), float(G3(p)), round(cluster_profile(3, float(p), DE, DI)[0][-1], 12)])
ok = [G3(F(11, 100)) < F(3, 2) < G3(F(3, 25)), G3(F(11, 50)) > F(19, 10)]
# p* with G3(p*) = 3/2 by exact bisection on rationals
lo, hi = F(11, 100), F(3, 25)
for _ in range(60):
    mid = (lo + hi) / 2
    lo, hi = (mid, hi) if G3(mid) < F(3, 2) else (lo, mid)
pstar = float(lo)
# baseline total cascade sizes Gamma = Gamma_C / (1 - Gamma_C rho D_I)
tot = lambda GC, rho: GC / (1 - GC * rho * DI)
rows += [['Gamma system 1 (n=1, rho=17/57)', str(tot(F(1), F(17, 57))), float(tot(F(1), F(17, 57))), ''],
         ['Gamma system 2 (n=3, G3(p*)=3/2, rho=5/38)', str(tot(F(3, 2), F(5, 38))), float(tot(F(3, 2), F(5, 38))), ''],
         ['Gamma system 3 (n=2, p=1-1/sqrt2, Gamma_C=3/2, rho=5/38)', str(tot(F(3, 2), F(5, 38))), float(tot(F(3, 2), F(5, 38))), round(cluster_profile(2, 1 - 1 / sqrt(2), DE, DI)[0][-1], 12)],
         ['Y, system 1: rho_Y D_I', str(2 * F(17, 57) * DI), float(2 * F(17, 57) * DI), ''],
         ['Y, system 2: lower bound G3(11/50) * 10/19', str(G3(F(11, 50)) * F(10, 19)), float(G3(F(11, 50)) * F(10, 19)), round(cluster_profile(3, 2 * pstar, DE, DI)[0][-1] * 10 / 19, 12)],
         ['Y, system 3: (2 sqrt2 - 1) * 10/19', '', (2 * sqrt(2) - 1) * 10 / 19, round(cluster_profile(2, 2 * (1 - 1 / sqrt(2)), DE, DI)[0][-1] * 10 / 19, 12)]]
rows.append(['all inequalities of the construction hold', str(all(ok) and 2 * F(17, 57) * DI > 1 and G3(F(11, 50)) * F(10, 19) > 1 and (2 * sqrt(2) - 1) * 10 / 19 < 1), '', ''])
write_csv('lemma_A73_exact.csv', ['quantity', 'exact', 'decimal', 'chain_cross_check'], rows)
