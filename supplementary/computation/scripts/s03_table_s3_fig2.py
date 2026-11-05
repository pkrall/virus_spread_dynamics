"""Table S3 and Fig. 2: replacement of X (at its fixed point) by Y with tau_Y = 1.5 or 2 tau_X in single TRTUs."""
import numpy as np
from _common import *
from trtu_model.trtu import fixed_point, simulate, replacement_times
rows = []; fig, axs = plt.subplots(1, 5, figsize=(10, 2.4), sharey=False)
Ts = 50
for k, (tau, K, m, DE, DI) in SCEN.items():
    row = [k, round(fixed_point(tau, K, m, DE, DI)[0], 4)]
    for f in (1.5, 2.0):
        r = simulate(tau, K, m, DE, DI, 4000, tauY=f * tau, TsY=Ts, start='fp')
        X0 = fixed_point(tau, K, m, DE, DI)[0]; Y0 = fixed_point(f * tau, K, m, DE, DI)[0]
        rt = replacement_times(r['NX'], r['NY'], X0, Ts)
        tot = (r['NX'] + r['NY']); imax = Ts + 1 + int(np.argmax(tot[Ts + 1:]))
        post_min = float(tot[imax:].min())
        row += [round(Y0, 4), rt['counted_1_95'], rt['elapsed_to_1'], round(float(tot[imax]), 4), round(post_min, 4)]
        if f == 2.0:
            ax = axs[k - 1]; t = np.arange(len(tot)) - Ts
            ax.plot(t[:350], r['NX'][:350], lw=1, label='X'); ax.plot(t[:350], r['NY'][:350], lw=1, label='Y'); ax.plot(t[:350], tot[:350], lw=0.8, ls='--', color='0.4', label='X+Y')
            ax.set_title(f'scenario {k}', loc='left'); ax.set_xlabel('steps after seeding of Y')
    rows.append(row)
axs[0].set_ylabel(r'new exposures'); axs[0].legend(frameon=False)
savefig(fig, 'fig2_strain_replacement_trtu.png')
hdr = ['scenario', 'p_N,X*']
for f in ('1.5', '2.0'):
    hdr += [f'p_N,Y*({f})', f'counted_1_to_95%({f})', f'elapsed_to_X<=1%({f})', f'max_total({f})', f'post-peak_min_total({f})']
write_csv('table_S3.csv', hdr, rows)
