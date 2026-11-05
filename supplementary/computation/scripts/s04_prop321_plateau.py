"""Proposition 3.2.1, numerical illustration and Fig. 3: three scenario-3 units seeded 10 steps apart."""
import numpy as np
from _common import *
from trtu_model.trtu import simulate
T = 400
I = simulate(0.12, 10, 0.035, 3, 4, T + 50)['IX']           # unit seeded at step 0
w, off = (0.4, 0.29, 0.31), (0, 10, 20)
comp = [np.concatenate([np.zeros(o), wi * I[:T - o]]) for wi, o in zip(w, off)]
y = sum(comp); mx = y.max(); im = int(np.argmax(y))
lo = im
while y[lo - 1] >= 0.95 * mx: lo -= 1
hi = im
while y[hi + 1] >= 0.95 * mx: hi += 1
def cross(level, a, b, up):
    for t in range(a, b):
        u, v = y[t], y[t + 1]
        if (u < level <= v) if up else (u >= level > v): return t + (level - u) / (v - u)
L_interp = cross(0.95 * mx, hi, hi + 3, False) - cross(0.95 * mx, lo - 3, lo, True)
tmin = hi + int(np.argmin(y[hi:hi + 60])); tmax2 = tmin + int(np.argmax(y[tmin:tmin + 80]))
tail = y[tmax2 + 28:]
rows = [['steps within 95% of maximum (integer)', hi - lo + 1], ['length with linear interpolation', round(L_interp, 2)],
        ['minimum after plateau (% of max)', round(100 * y[tmin] / mx, 1)], ['steps from last step >= 95% to minimum', tmin - hi],
        ['second maximum (% of max)', round(100 * y[tmax2] / mx, 1)], ['steps from minimum to second maximum', tmax2 - tmin],
        ['range from 28 steps after second maximum (% of max)', f'{100 * tail.min() / mx:.1f}-{100 * tail.max() / mx:.1f}'],
        ['limit (% of max)', round(100 * y[-1] / mx, 1)]]
write_csv('prop321_plateau_example.csv', ['quantity', 'value'], rows)
fig, (a1, a2) = plt.subplots(2, 1, figsize=(6.0, 4.6), sharex=True, gridspec_kw=dict(height_ratios=[1.1, 1]))
x = np.arange(T)
for c, wi, o in zip(comp, w, off): a1.plot(x, c, lw=1, label=f'unit seeded at step {o}, weight {wi}')
a1.plot(x, y, color='k', lw=1.6, label=r'aggregate $p_{I,\widehat{G},w}(t)$'); a1.set_ylabel('infectious mass'); a1.legend(frameon=False)
a2.plot(x, 100 * y / mx, color='k', lw=1.6); a2.axhspan(95, 100, color='0.85', lw=0); a2.axvspan(lo, hi, color='C0', alpha=0.12, lw=0)
a2.set_ylim(0, 105); a2.set_xlim(-12, 230); a2.set_ylabel('% of maximum'); a2.set_xlabel('steps after seeding of the first unit')
savefig(fig, 'fig3_plateau_trta.png')
