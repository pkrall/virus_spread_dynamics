"""Exact rational arithmetic for the three-system example in the proof of Proposition 3.3.7(a)
(D_E = D_I = 1), with a floating-point cross-check against the within-cluster chain."""
from fractions import Fraction as F
from _common import *
from trtu_model.cscp import cluster_profile
GXC = 1 + F(2, 5) + 2 * F(1, 25) - 2 * F(1, 125)        # n = 3, in-cluster probability 1/5
GYC3 = 1 + F(4, 5) + 2 * F(4, 25) - 2 * F(8, 125)       # n = 3, in-cluster probability 2/5
GYC2 = 1 + F(116, 125)                                  # n = 2, in-cluster probability 116/125
r = F(3, 8); rhoX = r / GXC; G = GXC / (1 - r)
rows = [['Gamma_X,C (n=3, 1/5)', str(GXC), float(GXC), cluster_profile(3, 0.2, 1, 1)[0][-1]],
        ['Gamma_Y,C (n=3, 2/5)', str(GYC3), float(GYC3), cluster_profile(3, 0.4, 1, 1)[0][-1]],
        ['Gamma_Y,C (n=2, 116/125)', str(GYC2), float(GYC2), cluster_profile(2, 116 / 125, 1, 1)[0][-1]],
        ['rho_X', str(rhoX), float(rhoX), ''], ['Gamma_X total', str(G), float(G), ''],
        ['R_Y system 2', str(GYC3 * 2 * rhoX), float(GYC3 * 2 * rhoX), ''], ['R_Y system 3', str(GYC2 * 2 * rhoX), float(GYC2 * 2 * rhoX), ''],
        ['R_Y system 1 (mass action)', str(2 * (G - 1) / G), float(2 * (G - 1) / G), '']]
write_csv('prop337a_exact.csv', ['quantity', 'exact', 'decimal', 'chain_cross_check'], rows)
