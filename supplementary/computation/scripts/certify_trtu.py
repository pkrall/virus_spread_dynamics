"""Computer-assisted proofs for Theorem 3.1.3: Table 1 (all five scenarios) and Tables S1/S2.

Writes to outputs/certificates/.  Standard library only (decimal); 80 significant digits,
outward rounding.  Run:  python scripts/certify_trtu.py
"""
import csv, os, sys, time
from decimal import Decimal as D
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from trtu_model.certify import (Unit, refine_root, naive_trajectory, dist_inf, contraction,
                                tail_range, certify_extrema)
from trtu_model.interval import I, fmt_dn, fmt_up
from trtu_model.trtu import fixed_point as fp_float

OUT = os.path.join(os.path.dirname(__file__), '..', 'outputs', 'certificates')
os.makedirs(OUT, exist_ok=True)
SCEN = {1: ('0.05', 8, '0.025', 3, 4), 2: ('0.05', 20, '0.035', 3, 4), 3: ('0.12', 10, '0.035', 3, 4),
        4: ('0.12', 20, '0.035', 3, 4), 5: ('0.12', 20, '0.02', 4, 5)}
TABLE1 = {1: (0.0071, 0.0127, 0.0058, 0.0076), 2: (0.0248, 0.0786, 0.0155, 0.0290), 3: (0.0266, 0.0932, 0.0164, 0.0312),
          4: (0.0311, 0.1473, 0.0216, 0.0347), 5: (0.0183, 0.1288, 0.0112, 0.0212)}
TABLE_S1 = {5: ((0.0458, 8), (0.0278, 23), (0.0321, 40)), 10: ((0.0455, 10), (0.0281, 25), (0.0321, 43)), 20: ((0.0397, 14), (0.0291, 32), (0.0317, 48))}
TABLE_S2 = {5: ((0.0171, 11), (0.0307, 37), (0.0253, 61)), 10: ((0.0178, 15), (0.0304, 40), (0.0254, 64)), 20: ((0.0197, 20), (0.0294, 48), (0.0257, 72))}
log = []
def say(*a):
    s = ' '.join(str(x) for x in a); print(s); log.append(s)

def r4(x):   # value rounded to 4 decimals as in the tables
    return float(round(float(x.mid()), 4))

# ------------------------------------------------------------ fixed points and contraction boxes
def float_power(tau, K, m, DE, DI):
    """Smallest P with |J^P|_inf < 0.9 for the floating-point Jacobian at the fixed point (a candidate only)."""
    import numpy as np
    q = 1 - m; d = DE + 1; L = DE + DI
    N, S, Ip = fp_float(tau, K, m, DE, DI); uu = 1 - tau * Ip; n = L + 1
    J = np.zeros((n, n)); J[0, 0] = q * uu ** K; J[1, 0] = 1 - uu ** K
    for j in range(d, L + 1):
        g = S * K * uu ** (K - 1) * tau * q ** j; J[0, j] = -q * g; J[1, j] = g
    for i in range(1, L): J[i + 1, i] = 1
    M = np.eye(n)
    for P in range(1, 2000):
        M = J @ M
        if np.abs(M).sum(1).max() < 0.9: return P
    return 24

def fixed_point_and_box(k):
    tau, K, m, DE, DI = SCEN[k]
    u = Unit(tau, K, m, DE, DI)
    g = fp_float(float(tau), K, float(m), DE, DI)[0]
    lo, hi = refine_root(u, g, 40)
    fp = u.fixed_point(lo, hi)
    assert fp['ok'], ('fixed point bracket failed', k)
    Pf = float_power(float(tau), K, float(m), DE, DI)
    for r in (D('5e-5'), D('1e-5'), D('1e-6'), D('1e-8')):
        for P in dict.fromkeys((Pf, 24, 48, 96)):
            c = contraction(u, fp['X'], r, P)
            if c['ok']:
                return u, fp, c, r, P
    raise RuntimeError(f'no contraction certificate for scenario {k}')

def run_until_box(step_fn, X1, c, r, t_min, t_max=20000):
    """Iterate X_{t+1} = step_fn(t, X_t) from X_1; stop at the first t >= t_min with |X_t - c| <= r."""
    Xs = [None, X1]; N = [None]
    t = 1
    while True:
        if t >= t_min and dist_inf(Xs[t], c) <= r:
            return Xs, N, t
        X = step_fn(t, Xs[t]); N.append(X[1]); Xs.append(X); t += 1
        if t > t_max: raise RuntimeError('trajectory did not enter the box')

def write_enclosures(name, N, Xs, T):
    with open(os.path.join(OUT, f'enclosures_{name}.csv'), 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['step', 'p_N_lower', 'p_N_upper', 'p_S_lower', 'p_S_upper'])
        for t in range(1 if N[0] is None else 0, len(N)):
            S = Xs[t][0] if t >= 1 else I(1)
            w.writerow([t, fmt_dn(N[t].lo), fmt_up(N[t].hi), fmt_dn(S.lo), fmt_up(S.hi)])

def write_tube(name, c):
    with open(os.path.join(OUT, f'tube_{name}.csv'), 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['box', 'coordinate', 'lower', 'upper'])
        labels = ['p_S(t)'] + [f'p_N(t-{j})' for j in range(1, len(c['B']))]
        for k, box in enumerate(c['tube']):
            for lab, x in zip(labels, box):
                w.writerow([k, lab, fmt_dn(x.lo), fmt_up(x.hi)])

t0 = time.time()
boxes = {}
fp_rows, con_rows = [], []
for k in SCEN:
    u, fp, c, r, P = fixed_point_and_box(k)
    boxes[k] = (u, fp, c, r, P)
    tr = tail_range(c['tube'])
    fp_rows.append([k, *SCEN[k], fmt_dn(fp['N'].lo, 40), fmt_up(fp['N'].hi, 40), fmt_dn(fp['H_lo'].lo, 12), fmt_up(fp['H_hi'].hi, 12), fp['ok']])
    con_rows.append([k, str(r), P, fmt_up(c['L'], 12), c['ok'], fmt_dn(tr.lo, 12), fmt_up(tr.hi, 12)])
    write_tube(f'scenario{k}', c)
    say(f'scenario {k}: fixed point p_N* in [{fmt_dn(fp["N"].lo, 20)}, {fmt_up(fp["N"].hi, 20)}]; box radius {r}, P = {P}, |DF^P| <= {fmt_up(c["L"], 8)} ({time.time() - t0:.0f}s)')
with open(os.path.join(OUT, 'fixed_points.csv'), 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['scenario', 'tau', 'K', 'm', 'D_E', 'D_I', 'p_N*_lower', 'p_N*_upper', 'H(lower)_lower_bound', 'H(upper)_upper_bound', 'sign_change'])
    w.writerows(fp_rows)
with open(os.path.join(OUT, 'contraction.csv'), 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['scenario', 'box_radius_r', 'power_P', 'upper_bound_|DF^P|_inf_on_B', 'F^P(B)_in_B_and_contracting', 'tail_p_N_lower', 'tail_p_N_upper'])
    w.writerows(con_rows)

# ------------------------------------------------------------ the bracket printed in Appendix A.1.3
u3 = boxes[3][0]
fpA = u3.fixed_point('0.02664168989885590', '0.02664168989885591')
say(f'Appendix A.1.3 bracket [0.02664168989885590, 0.02664168989885591]: H(lower) in {fpA["H_lo"]}, H(upper) in {fpA["H_hi"]}, sign change {fpA["ok"]}')

# ------------------------------------------------------------ Table 1
rows = []
for k in SCEN:
    u, fp, c, r, P = boxes[k]
    sig = I(D('0.01')) * I(u.m)
    X1 = [u.q * (I(1) - sig) + I(u.m), sig] + [I(0)] * (u.L - 1)
    Xs, N, T = run_until_box(lambda t, X: u.F(X), X1, c['center'], r, t_min=1)
    N[0] = sig
    tr = tail_range(c['tube'])
    ext = certify_extrema(N, list(range(1, T)), ('max', 'min', 'max'), tr)
    write_enclosures(f'table1_scenario{k}', N, Xs, T)
    exp = TABLE1[k]
    fp_ok = r4(fp['N']) == exp[0]
    rows.append([k, 'p_N*', '', fmt_dn(fp['N'].lo, 20), fmt_up(fp['N'].hi, 20), r4(fp['N']), exp[0], fp['ok'] and fp_ok])
    for e, val in zip(ext, exp[1:]):
        rows.append([k, {'max': 'maximum', 'min': 'lowest minimum after maximum'}[e['kind']] if e is ext[0] or e['kind'] == 'min' else 'highest maximum after that minimum',
                     e['step'], fmt_dn(e['enclosure'].lo, 20), fmt_up(e['enclosure'].hi, 20), r4(e['enclosure']), val, e['ok'] and r4(e['enclosure']) == val])
    say(f'Table 1 scenario {k}: box entered at step {T}; extrema at steps {[e["step"] for e in ext]}, values {[r4(e["enclosure"]) for e in ext]}, certified {all(e["ok"] for e in ext)}')
    if k == 3:
        say(f'  scenario 3: |X_203 - c|_inf <= {fmt_up(dist_inf(Xs[203], c["center"]), 12)} (box radius {r})' if T <= 203 else '')
with open(os.path.join(OUT, 'table1_certificate.csv'), 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['scenario', 'quantity', 'step', 'enclosure_lower', 'enclosure_upper', 'rounded', 'value_in_Table_1', 'certified_and_matches'])
    w.writerows(rows)

# ------------------------------------------------------------ Tables S1 and S2 (linear ramps of K)
rows = []
for table, K0, K1, s0, s1, kinds, ref in (('S1', 10, 20, 3, 4, ('max', 'min', 'max'), TABLE_S1), ('S2', 20, 10, 4, 3, ('min', 'max', 'min'), TABLE_S2)):
    u0, fp0 = boxes[s0][0], boxes[s0][1]
    u1, fp1, c1, r1, P1 = boxes[s1]
    for n in (5, 10, 20):
        def step(t, X, n=n):
            K = D(K0) + D(K1 - K0) * D(min(t, n)) / D(n)
            return u1.F(X, K)
        Xs, N, T = run_until_box(step, list(fp0['X']), c1['center'], r1, t_min=n + 1)
        tr = tail_range(c1['tube'])
        ext = certify_extrema(N, list(range(1, T)), kinds, tr)
        write_enclosures(f'table{table}_ramp{n}', N, Xs, T)
        for e, (val, st) in zip(ext, ref[n]):
            rows.append([table, n, e['kind'], e['step'], st, fmt_dn(e['enclosure'].lo, 20), fmt_up(e['enclosure'].hi, 20), r4(e['enclosure']), val,
                         e['ok'] and r4(e['enclosure']) == val and e['step'] == st])
        say(f'Table {table}, {n}-step ramp: box entered at step {T}; extrema {[(e["kind"], e["step"], r4(e["enclosure"])) for e in ext]}, certified {all(e["ok"] for e in ext)}')
with open(os.path.join(OUT, 'tables_S1_S2_certificate.csv'), 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['table', 'ramp_steps', 'kind', 'step', 'step_in_table', 'enclosure_lower', 'enclosure_upper', 'rounded', 'value_in_table', 'certified_and_matches'])
    w.writerows(rows)
say(f'total time {time.time() - t0:.0f}s')
open(os.path.join(OUT, 'summary.txt'), 'w').write('\n'.join(log) + '\n')
