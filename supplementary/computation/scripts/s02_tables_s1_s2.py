"""Tables S1, S2 and Figs. S2, S3: linear ramps of K between scenarios 3 and 4, starting at the fixed point."""
import numpy as np
from _common import *
from trtu_model.trtu import fixed_point, ramp, successive_extrema
rows = []
for table, K0, K1, kinds, fname, title in (('S1', 10, 20, ('max', 'min', 'max'), 'figS2_K_ramp_up.png', 'K ramped from 10 to 20'),
                                           ('S2', 20, 10, ('min', 'max', 'min'), 'figS3_K_ramp_down.png', 'K ramped from 20 to 10')):
    fig, ax = plt.subplots(figsize=(3.0, 3.1))
    for n in (5, 10, 20):
        N = ramp(0.12, K0, K1, n, 0.035, 3, 4, T=600)
        ext = successive_extrema(N, kinds, start=0)
        rows.append([table, n, round(fixed_point(0.12, K1, 0.035, 3, 4)[0], 4)] + sum([[round(v, 4), i + 1] for v, i in ext], []))
        ax.plot(np.arange(1, 151), N[:150], lw=1.1, label=f'{n}-step ramp')
    ax.set_title(title, loc='left'); ax.set_xlabel('steps after start at the fixed point'); ax.set_ylabel(r'new exposures $p_N(t)$')
    ax.set_xlim(0, 150); ax.legend(frameon=False); savefig(fig, fname)
write_csv('tables_S1_S2.csv', ['table', 'ramp_steps', 'p_N*_final', 'extremum1', 'step1', 'extremum2', 'step2', 'extremum3', 'step3'], rows)
