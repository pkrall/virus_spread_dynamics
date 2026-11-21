"""Table 4 and the remark (scale of the approximation) after Theorem 5: the finite-population variant of Section 3.2,
solved numerically for one transmission unit of scenario 1 of Table 2 (tau_Y = 2 tau_X) coupled to the three CSCPs
of Table 3 (alpha = 0.1, population of the transmission unit normalized to N_U = 1).  Finite start of the remark after
Lemma A.3.1: X is seeded into a naive transmission unit at step 0, the seeded population is empty and no clusters
are affected (C_0 = D_0 = 0), Y is seeded at step 400, and the run ends at step 1100.

Variant rules (Section 3.2): (i) departures join the seeded population; departed non-susceptibles are counted in D;
(ii) arrivals are exposed/infectious with V with probability a_V(t) and immune with probability z(t); the
susceptible inflow is m(1 - a_X - a_Y - z), m z enters the removed compartment and m a_V/(1-m) is added to the
cohort of strain V entering the exposed chain; (iii) seedings are thinned by 1 - u(t-1), and u increases by
r n times the number of newly seeded clusters.  r = N_U/N_C; r = 0 is the model.
Output (table4_theorem5_scale.csv), with and without infected arrivals: counted steps of the Y share from 1% to 95%,
elapsed steps from the introduction of Y until the exposure intensity of X falls to 1% of the X-only stationary
level of the model, elapsed steps from 90% to 1% of that level (empty if X is below 90% already when Y is
introduced), steps until the exposures of X in the transmission unit fall to 1% of their level before Y, the
exposure intensity of X at the introduction of Y in % of the stationary level, the non-susceptible fraction of the
seeded population at step 400 and its largest value over the run, and the largest deviation of the compartment
sum of the transmission unit from one."""
import numpy as np
from _common import *
from trtu_model.trtu import fixed_point
from trtu_model.cscp import cluster_profile, solve_rb, total_profile

tau, K, m, DE, DI = SCEN[1]; q, d, L = 1 - m, DE + 1, DE + DI; alpha = 0.1
rb2 = solve_rb(7, DE, DI, 4.0); rb3 = solve_rb(5, DE, DI, 3.0)
SYS = {'x1': (1, 0.0, 0.225), 'x2': (7, rb2, 0.0375), 'x3': (5, rb3, 0.7 / 12)}


def run(name, r, T=1100, TsY=400, imports=True):
    n, rb, rho = SYS[name]
    G, ib, g = cluster_profile(n, rb, DE, DI); H = len(g)
    e = np.array([sum(g[a - k] for k in range(DE) if a - k >= 0) for a in range(H)]); rem = G - e - ib
    S = np.ones(T); cX = np.zeros(T); cY = np.zeros(T); pNX = np.zeros(T); pNY = np.zeros(T)
    IXs = np.zeros(T); IYs = np.zeros(T); JX = np.zeros(T); JY = np.zeros(T); IX = np.zeros(T); IY = np.zeros(T)
    a = np.zeros((2, T)); z = np.zeros(T); u = np.zeros(T); D = np.zeros(T); R = np.zeros(T); err = 0.0
    sh = lambda arr, j, t: arr[t - j] if t - j >= 0 else 0.0
    conv = lambda prof, J, t: sum(prof[k] * J[t - k] for k in range(min(H, t + 1)))
    for t in range(T):
        if t > 0:
            S[t] = (S[t - 1] - pNX[t - 1] - pNY[t - 1]) * q + m * (1 - a[0, t - 1] - a[1, t - 1] - z[t - 1])
            D[t] = D[t - 1] + m * (1 - S[t - 1] + pNX[t - 1] + pNY[t - 1])
            R[t] = R[t - 1] * q + (sh(cX, L + 1, t) + sh(cY, L + 1, t)) * q ** (L + 1) + m * z[t - 1]
        IXs[t] = sum(q ** j * sh(cX, j, t) for j in range(d, L + 1)); IYs[t] = sum(q ** j * sh(cY, j, t) for j in range(d, L + 1))
        E = sum(q ** j * (sh(cX, j, t) + sh(cY, j, t)) for j in range(1, DE + 1))
        err = max(err, abs(S[t] + E + IXs[t] + IYs[t] + R[t] - 1))
        w = tau * IXs[t] + 2 * tau * IYs[t]; lam = 1 - (1 - w) ** K
        nx = lam * S[t] * tau * IXs[t] / w if w > 0 else 0.0
        ny = lam * S[t] * 2 * tau * IYs[t] / w if w > 0 else 0.0
        sx = 0.01 * m if t == 0 else 0.0; sy = 0.01 * m if t == TsY else 0.0
        if sx + sy > 0:
            nx = nx * (1 - (sx + sy) / S[t]) + sx; ny = ny * (1 - (sx + sy) / S[t]) + sy
        pNX[t] = nx; pNY[t] = ny
        up = u[t - 1] if t > 0 else 0.0
        crossX = rho * sum(ib[k] * JX[t - 1 - k] for k in range(min(H, t))); crossY = rho * sum(ib[k] * JY[t - 1 - k] for k in range(min(H, t)))
        JX[t] = (1 - up) * (alpha * IXs[t] + crossX); JY[t] = (1 - up) * (alpha * IYs[t] + crossY)
        u[t] = up + r * n * (JX[t] + JY[t])
        IX[t] = conv(g, JX, t); IY[t] = conv(g, JY, t)
        old = JX[:max(0, t - H + 1)].sum() + JY[:max(0, t - H + 1)].sum()
        a[0, t] = r * conv(e + ib, JX, t) * imports; a[1, t] = r * conv(e + ib, JY, t) * imports
        z[t] = r * (conv(rem, JX, t) + conv(rem, JY, t) + G[-1] * old + D[t])
        cX[t] = nx + m * a[0, t] / q; cY[t] = ny + m * a[1, t] / q
    return dict(pNX=pNX, IX=IX, IY=IY, ns=(a.sum(0) + z) if imports else (r * 0 + z), err=err)


_, _, IX0 = fixed_point(tau, K, m, DE, DI)
Ts = 400


def measures(name, o):
    n, rb, rho = SYS[name]
    I0 = alpha * IX0 * total_profile(n, rb, rho, DE, DI, 3000)['Gamma_tot']
    IX, IY = o['IX'], o['IY']
    x1 = next(i for i in range(Ts + 1, len(IX)) if IX[i] <= 0.01 * I0) - Ts
    x90 = next(i for i in range(Ts + 1, len(IX)) if IX[i] <= 0.9 * I0) - Ts
    sh_ = IY / np.where(IX + IY > 0, IX + IY, 1)
    t1 = next(i for i in range(Ts + 1, len(IX)) if sh_[i] >= 0.01); t95 = next(i for i in range(Ts + 1, len(IX)) if sh_[i] >= 0.95)
    src = next(i for i in range(Ts + 1, len(IX)) if o['pNX'][i] <= 0.01 * o['pNX'][Ts - 1]) - Ts
    x90_1 = x1 - x90 if IX[Ts] > 0.9 * I0 else ''
    return t95 - t1 + 1, x1, x90_1, src, IX[Ts] / I0


rows = []
for name, (n, rb, rho) in SYS.items():
    for label, r in (('0', 0.0), ('2e-5', 2e-5), ('2e-4', 2e-4), ('1e-3', 1e-3), ('1/(100n)', 1 / (100 * n))):
        o = run(name, r); on = run(name, r, imports=False)
        ys, x1, x90_1, src, lvl = measures(name, o); ysn, x1n, x90_1n, _, _ = measures(name, on)
        rows.append([name, label, f'{r:.6g}', ys, x1, x90_1, ysn, x1n, x90_1n, src, round(100 * lvl, 1),
                     round(float(o['ns'][Ts]), 4), round(float(o['ns'].max()), 4), f"{o['err']:.1e}"])
write_csv('table4_theorem5_scale.csv', ['system', 'r_label', 'r', 'Y_share_1_to_95%_counted', 'X_to_1%_elapsed', 'X_90_to_1%_elapsed',
                                        'without_infected_arrivals_Y_share', 'without_infected_arrivals_X_to_1%', 'without_infected_arrivals_X_90_to_1%',
                                        'source_steps_X_to_1%', 'X_at_introduction_of_Y_%_of_stationary_level',
                                        'non-susceptible_fraction_of_seeded_population_at_step_400', 'max_non-susceptible_fraction_steps_0_to_1100',
                                        'max_|compartment_sum-1|'], rows)
