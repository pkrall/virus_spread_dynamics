"""Remarks after Lemma A.1.2 (return transmission) and after Lemma A.3.1 (finite-start accounting).

(1) Return transmission, an experiment outside the model: exposed arrivals into the transmission unit of each strain V
    proportional to its exposure intensity in the seeded population, sigma_V(t) = phi * I_V(t-1),
    phi = 0, 0.001, 0.01, 0.1.  Reports the exposed share of the turnover inflow at the X-only steady state and
    the steps until X falls to 1% of its previous level in the transmission unit and in the seeded population.
(2) Finite start: X seeded into a naive transmission unit at step 0 with an empty seeded population, Y seeded at
    step 400.  Reports clusters seeded before Y, clusters and exposures per individual in the transmission units over the window.
Transmission unit: scenario 1 of Table 2, tau_Y = 2 tau_X; seeded populations: the three CSCPs of Table 3; alpha/N_U = 0.1."""
import numpy as np
from _common import *
from trtu_model.trtu import fixed_point, simulate
from trtu_model.cscp import solve_rb, total_profile

tau, K, m, DE, DI = SCEN[1]; alpha = 0.1; G = 10
q, d, L = 1 - m, DE + 1, DE + DI
N0, S0, IX0 = fixed_point(tau, K, m, DE, DI)
rb2 = solve_rb(7, DE, DI, 4.0); rb3 = solve_rb(5, DE, DI, 3.0)
sysp = {'x1': (1, 0.0, 9 / 40), 'x2': (7, rb2, (G - 4) / (G * 4 * DI)), 'x3': (5, rb3, (G - 3) / (G * 3 * DI))}
prof = {k: total_profile(n, rb, rho, DE, DI, A=4000) for k, (n, rb, rho) in sysp.items()}


def run(phi, h, T1=6000, T2=1500):
    n = T1 + T2; hh = h[:1500]
    ntx = np.zeros(n); nty = np.zeros(n); s = np.ones(n); kx = np.zeros(n); ky = np.zeros(n)
    Ix = np.zeros(n); Iy = np.zeros(n); nx = np.zeros(n)
    ntx[:L + 2] = N0; s[:L + 2] = S0
    for t in range(L + 2, n):
        s[t] = (s[t - 1] - ntx[t - 1] - nty[t - 1]) * q + m
        ix = sum(q ** j * ntx[t - j] for j in range(d, L + 1)); iy = sum(q ** j * nty[t - j] for j in range(d, L + 1))
        kx[t] = alpha * ix; ky[t] = alpha * iy
        a0 = max(0, t - len(hh) + 1)
        Ix[t] = np.dot(kx[a0:t + 1][::-1], hh[:t - a0 + 1]); Iy[t] = np.dot(ky[a0:t + 1][::-1], hh[:t - a0 + 1])
        w = ix * tau + iy * 2 * tau; lam = 1 - (1 - w) ** K
        ex = lam * s[t] * ix * tau / w if ix > 0 else 0.0; ey = lam * s[t] * iy * 2 * tau / w if iy > 0 else 0.0
        sx = phi * Ix[t - 1]; sy = phi * Iy[t - 1] + (0.01 * m if t == T1 else 0.0)
        nx[t] = ex; ntx[t] = ex * (1 - (sx + sy) / s[t]) + sx; nty[t] = ey * (1 - (sx + sy) / s[t]) + sy
    ref = nx[T1 - 1]; refI = Ix[T1 - 1]
    t_src = next((i for i in range(1, T2) if nx[T1 + i] <= 0.01 * ref), None)
    t_snk = next((i for i in range(1, T2) if Ix[T1 + i] <= 0.01 * refI), None)
    return phi * refI / m, t_src, t_snk


rows = []
for k in sysp:
    for phi in (0, 1e-3, 1e-2, 1e-1):
        share, ts, tk = run(phi, prof[k]['g'])
        rows.append([k, phi, round(100 * share, 2), ts, tk])
write_csv('remark_return_transmission.csv', ['system', 'phi', 'exposed_share_of_turnover_inflow_%', 'source_steps_X_to_1%', 'seeded_population_steps_X_to_1%'], rows)

TY = 400; r = simulate(tau, K, m, DE, DI, TY + 1200, tauY=2 * tau, TsY=TY, start='naive', TsX=0)
kX = alpha * r['IX']; kY = alpha * r['IY']
rows = []
for k in sysp:
    P = prof[k]; h = P['g']; B = P['B']; I0 = alpha * IX0 * P['Gamma_tot']
    IX = np.convolve(h, kX)[:len(kX)]; IY = np.convolve(h, kY)[:len(kX)]
    x1 = next(i for i in range(TY + 1, len(IX)) if IX[i] <= 0.01 * I0) - TY
    J = np.convolve(B, kX + kY)[:len(kX)]          # newly seeded clusters per step (direct and cross-cluster)
    rows.append([k, round(IX[TY] / I0, 4), x1, f'0-{TY + x1}', round(J[:TY + 1].sum(), 2), round(J[:TY + x1 + 1].sum(), 2), round((IX + IY)[:TY + x1 + 1].sum(), 1)])
write_csv('remark_finite_start.csv', ['system', 'I_X_at_step_400/stationary', 'steps_after_400_until_X<=1%', 'accounting_window', 'clusters_seeded_through_step_400', 'clusters_seeded_over_window', 'seeded_population_exposures_over_window'], rows)
