"""Numbers in the remarks on rates and shares (Appendix A):
- remark on the rate of elimination after Lemma A.1.2: source bound (tau_X/tau_Y)^(1/L) and the numerical limit of
  the growth factor p_N,X(t+1)/p_N,X(t), scenario 1, tau_Y = 2 tau_X;
- Corollary A.4.1(f): decay factors lambda of the response kernels of the three CSCPs of Table 3, i.e. the root
  in (0, 1) of sum_a beta(a) x^(-a) = 1, beta the cross-cluster seeding kernel;
- Lemma A.2.2(a): the bound 2A/B for Gamma = 10, n = 7, D_E = 3, D_I = 4, and the limiting ratio of the
  exposure-reduction times, ln(lambda_mass action) / ln(lambda_clustered) as the in-cluster probability tends to 1;
- remark after Lemma A.3.1: long-run shares s(Gamma-1)/(1+s Gamma) and s Gamma/(1+s Gamma) for s = (alpha/N_U) A_I
  and s = 1;
- remark (finite populations) after Lemma A.3.1: over the replacement windows of Table 3 (stationary prehistory,
  Y seeded at step 50), the share of all exposures generated within the CSCP, exposures per individual in the transmission units,
  and newly seeded clusters per individual in the transmission units (direct and cross-cluster)."""
import numpy as np
from _common import *
from trtu_model.trtu import fixed_point, simulate, replacement_times
from trtu_model.cscp import solve_rb, total_profile, cluster_profile, response

tau, K, m, DE, DI = SCEN[1]; L = DE + DI; alpha = 0.1; G = 10
rows = []
# source rate
src_bound = 0.5 ** (1 / L)
rr = simulate(tau, K, m, DE, DI, 2000, tauY=2 * tau, TsY=50, start='fp')
gf = rr['NX'][1501] / rr['NX'][1500]
rows += [['source rate bound (tau_X/tau_Y)^(1/L)', round(src_bound, 4)], ['growth factor p_N,X(t+1)/p_N,X(t) at t = 1500', round(gf, 4)]]


def decay_factor(n, rb, rho):
    _, ib, _ = cluster_profile(n, rb, DE, DI)
    beta = rho * ib                                 # cross-cluster seedings at cascade age a + 1
    F = lambda x: sum(b * x ** -(a + 1) for a, b in enumerate(beta)) - 1
    lo, hi = 1e-9, 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if F(mid) > 0 else (lo, mid)
    return 0.5 * (lo + hi)


rb2 = solve_rb(7, DE, DI, 4.0); rb3 = solve_rb(5, DE, DI, 3.0)
sysp = {'x1': (1, 0.0, 0.225), 'x2': (7, rb2, 0.0375), 'x3': (5, rb3, 0.7 / 12)}
lam = {k: decay_factor(*v) for k, v in sysp.items()}
for k in sysp: rows.append([f'decay factor lambda_X, CSCP {k}', round(lam[k], 4)])
# Lemma A.2.2(a)
A = np.log(G / (G - 1)); Bv = np.log(G / (G - 7))
rows.append(['2A/B, Gamma = 10, n = 7', round(2 * A / Bv, 4)])
lam_lim = decay_factor(7, 1.0, (G - 7) / (G * 7 * DI))
rows.append(['limiting ratio ln(lambda_x1)/ln(lambda_n=7, in-cluster probability 1)', round(np.log(lam['x1']) / np.log(lam_lim), 4)])
for rb in (0.9, 0.99, 0.999):
    GC = cluster_profile(7, rb, DE, DI)[0][-1]
    rows.append([f'  ratio at in-cluster probability {rb}', round(np.log(lam['x1']) / np.log(decay_factor(7, rb, (G - GC) / (G * GC * DI))), 4)])
# shares
_, _, IX0 = fixed_point(tau, K, m, DE, DI); NX0 = fixed_point(tau, K, m, DE, DI)[0]
AI = sum((1 - m) ** j for j in range(DE + 1, L + 1)); s = alpha * AI
rows.append(['A_I, scenario 1', round(AI, 4)]); rows.append(['s = (alpha/N_U) A_I', round(s, 4)])
for sv in (s, 1.0):
    rows.append([f'share generated within CSCP, s = {sv:.4f}', round(sv * (G - 1) / (1 + sv * G), 4)])
    rows.append([f'share occurring in CSCP, s = {sv:.4f}', round(sv * G / (1 + sv * G), 4)])
write_csv('rates_and_shares.csv', ['quantity', 'value'], rows)
# replacement windows
Ts = 50; T = 2500; pre = 3000
rr = simulate(tau, K, m, DE, DI, T, tauY=2 * tau, TsY=Ts, start='fp')
rows = []
for k, (n, rb, rho) in sysp.items():
    P = total_profile(n, rb, rho, DE, DI, 3000)
    kx = np.concatenate([np.full(pre, alpha * IX0), alpha * rr['IX']]); ky = np.concatenate([np.zeros(pre), alpha * rr['IY']])
    IXc = response(kx, P['g'])[pre:]; IYc = response(ky, P['g'])[pre:]
    J = np.convolve(P['B'], kx + ky)[:len(kx)][pre:]
    w = replacement_times(IXc, IYc, alpha * IX0 * P['Gamma_tot'], Ts)['elapsed_to_1']
    win = slice(Ts + 1, Ts + w + 1)
    src = (rr['NX'] + rr['NY'])[win].sum(); csc = (IXc + IYc)[win].sum(); seeds = alpha * (rr['IX'] + rr['IY'])[win].sum()
    rows.append([k, w, round(csc, 1), round((csc - seeds) / (src + csc), 3), round(J[win].sum(), 2)])
write_csv('replacement_window_accounting.csv', ['system', 'window_steps', 'CSCP_exposures_per_source_individual', 'share_of_all_exposures_generated_within_CSCP', 'clusters_newly_seeded_per_source_individual'], rows)
