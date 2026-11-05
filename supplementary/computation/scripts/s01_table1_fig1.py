"""Table 1 and Fig. 1: single-strain TRTUs seeded with 0.01 m at step 0 (floating point)."""
import numpy as np
from _common import *
from trtu_model.trtu import fixed_point, simulate, successive_extrema
rows = []; fig, ax = plt.subplots(figsize=(5.5, 3.2))
for k, (tau, K, m, DE, DI) in SCEN.items():
    N = simulate(tau, K, m, DE, DI, 3000)['NX']
    (mx, imx), (mn, imn), (mx2, imx2) = successive_extrema(N, ('max', 'min', 'max'), start=1)
    Ns = fixed_point(tau, K, m, DE, DI)[0]
    rows.append([k, tau, K, m, DE, DI, round(Ns, 4), round(mx, 4), imx, round(mn, 4), imn, round(mx2, 4), imx2])
    ax.plot(np.arange(200), N[:200], lw=1, label=f'{k}: τ={tau}, K={K}, D_E={DE}, D_I={DI}, m={m}')
write_csv('table1.csv', ['scenario', 'tau', 'K', 'm', 'D_E', 'D_I', 'p_N*', 'max', 'step_max', 'lowest_min_after_max', 'step_min', 'highest_max_after_min', 'step_max2'], rows)
ax.set_xlabel('step after seeding'); ax.set_ylabel(r'new exposures $p_N(t)$'); ax.legend(frameon=False)
savefig(fig, 'fig1_trtu_exposure_curves.png')
