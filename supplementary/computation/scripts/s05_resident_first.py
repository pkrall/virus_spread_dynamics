"""Numerical illustration after Proposition 3.2.1 and Fig. S4: plateau of Y while it replaces an established X
(two scenario-1 units, tau_Y = 2 tau_X, Y introduced 18 steps apart, weights 0.6 and 0.4)."""
import numpy as np
from _common import *
from trtu_model.trtu import simulate, fixed_point, aggregate
tau, K, m, DE, DI = SCEN[1]
r = simulate(tau, K, m, DE, DI, 900, tauY=2 * tau, TsY=0, start='fp')
w, off = (0.6, 0.4), (0, 18)
y = aggregate(r['IY'], w, off, 800, pre=0.0); x = aggregate(r['IX'], w, off, 800)
Xs = fixed_point(tau, K, m, DE, DI)[2]
mx = y.max(); im = int(np.argmax(y)); lo = im
while y[lo - 1] >= 0.95 * mx: lo -= 1
hi = im
while y[hi + 1] >= 0.95 * mx: hi += 1
tmin = hi + int(np.argmin(y[hi:hi + 60])); tmax2 = tmin + int(np.argmax(y[tmin:tmin + 80]))
t1 = int(np.argmax(x < 0.01 * Xs))
rows = [['steps of Y within 95% of its maximum', hi - lo + 1], ['window (steps after first introduction)', f'{lo}-{hi}'],
        ['X at window start (% of endemic level)', round(100 * x[lo] / Xs, 1)], ['X at window end (% of endemic level)', round(100 * x[hi] / Xs, 1)],
        ['steps from window end until X < 1% of endemic level', t1 - hi], ['Y minimum after window (% of max)', round(100 * y[tmin] / mx, 1)],
        ['steps from window end to minimum', tmin - hi], ['Y rebound maximum (% of max)', round(100 * y[tmax2] / mx, 1)],
        ['Y limit (% of max)', round(100 * y[-1] / mx, 1)]]
write_csv('prop321_resident_first_example.csv', ['quantity', 'value'], rows)
fig, ax = plt.subplots(figsize=(6.0, 3.0)); t = np.arange(len(y))
ax.plot(t, y, color='C3', lw=1.6, label='strain Y (aggregate)'); ax.plot(t, x, color='C0', lw=1.4, label='strain X (aggregate)')
ax.plot(t, x + y, color='0.4', lw=0.9, ls='--', label='X+Y'); ax.axvspan(lo, hi, color='C3', alpha=0.1, lw=0)
ax.set_xlim(-5, 260); ax.set_xlabel('steps after the first introduction of Y'); ax.set_ylabel('infectious mass'); ax.legend(frameon=False)
savefig(fig, 'figS4_plateau_during_replacement.png')
