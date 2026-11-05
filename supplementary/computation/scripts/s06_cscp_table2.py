"""Table 2, Figs. 4 and 5: one scenario-1 source unit coupled (alpha = 0.1) to three CSCPs with total cascade size 10.
x1: mass action (n = 1, rho = 0.225); x2: n = 7, Gamma_C = 4, rho = 0.0375; x3: n = 5, Gamma_C = 3, rho = 0.7/12."""
import numpy as np
from _common import *
from trtu_model.trtu import fixed_point, simulate, successive_extrema, replacement_times
from trtu_model.cscp import solve_rb, total_profile, response
tau, K, m, DE, DI = SCEN[1]; alpha = 0.1; A = 3000; T = 2500
r2 = solve_rb(7, DE, DI, 4.0); r3 = solve_rb(5, DE, DI, 3.0)
print('in-cluster probabilities: x2 %.11f, x3 %.11f' % (r2, r3))
prof = {'x1': total_profile(1, 0.0, 0.225, DE, DI, A), 'x2': total_profile(7, r2, 0.0375, DE, DI, A), 'x3': total_profile(5, r3, 0.7 / 12, DE, DI, A)}
rows = []
# single strain: naive source unit seeded at step 0
src = simulate(tau, K, m, DE, DI, T)
(mx, _), (mn, _), (mx2, _) = successive_extrema(src['NX'], ('max', 'min', 'max'))
NX0, _, IX0 = fixed_point(tau, K, m, DE, DI); NY0, _, IY0 = fixed_point(2 * tau, K, m, DE, DI)
# replacement: X at its fixed point, Y seeded at step Ts with tau_Y = 2 tau_X
Ts = 50; rr = simulate(tau, K, m, DE, DI, T, tauY=2 * tau, TsY=Ts, start='fp')
rt = replacement_times(rr['NX'], rr['NY'], NX0, Ts)
rows.append(['TRTU (p_N)', round(NX0, 4), round(mx, 4), round(mn, 4), round(mx2, 4), round(NY0, 4), rt['counted_1_95'], rt['elapsed_to_1'], rt['elapsed_90_1']])
fig4, ax4 = plt.subplots(figsize=(4.5, 3.0)); fig5, ax5 = plt.subplots(figsize=(4.5, 3.0))
pre = A
for name, P in prof.items():
    g = P['g']
    I = response(alpha * src['IX'], g)
    (a, _), (b, _), (c, _) = successive_extrema(I, ('max', 'min', 'max'))
    kx = np.concatenate([np.full(pre, alpha * IX0), alpha * rr['IX']]); ky = np.concatenate([np.zeros(pre), alpha * rr['IY']])
    IXc = response(kx, g)[pre:]; IYc = response(ky, g)[pre:]
    X0 = alpha * IX0 * P['Gamma_tot']
    rc = replacement_times(IXc, IYc, X0, Ts)
    rows.append([f'CSCP {name} (I)', round(I[-1], 4), round(a, 4), round(b, 4), round(c, 4), round(alpha * IY0 * P['Gamma_tot'], 4), rc['counted_1_95'], rc['elapsed_to_1'], rc['elapsed_90_1']])
    ax4.plot(I[:300], lw=1.1, label=name); ax5.plot(np.arange(len(IYc[:500])) - Ts, (IYc / (IXc + IYc))[:500], lw=1.1, label=name)
write_csv('table2.csv', ['system', 'steady_state', 'max', 'lowest_min_after_max', 'highest_max_after_min', 'Y_steady_state', 'counted_Y_share_1_to_95%', 'elapsed_to_X<=1%', 'elapsed_X_90%_to_1%'], rows)
write_csv('table2_profiles.csv', ['system', 'n', 'rb', 'rho', 'Gamma_C', 'Gamma_total'], [[k, n, rb, rho, round(prof[k]['GammaC'], 8), round(prof[k]['Gamma_tot'], 8)] for k, n, rb, rho in (('x1', 1, 0.0, 0.225), ('x2', 7, r2, 0.0375), ('x3', 5, r3, 0.7 / 12))])
ax4.set_xlabel('step after seeding of the source unit'); ax4.set_ylabel(r'exposure intensity $\mathcal{I}(t)$'); ax4.legend(frameon=False); savefig(fig4, 'fig4_cscp_fluctuations.png')
ax5.set_xlabel('steps after seeding of Y in the source'); ax5.set_ylabel('Y share of CSCP exposures'); ax5.legend(frameon=False); savefig(fig5, 'fig5_cscp_replacement.png')
