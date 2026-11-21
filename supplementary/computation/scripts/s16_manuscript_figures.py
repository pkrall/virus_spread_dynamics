"""Figures 2-5 and 7 of the article in their published layout (Fig. 1, the model schematic, is written by s18,
Fig. 6 by s04, the supplementary figures by s02 and s05).  Output names follow the article: fig2.png ... fig7.png."""
import numpy as np
from matplotlib.lines import Line2D
from _common import *
from trtu_model.trtu import fixed_point, simulate
from trtu_model.cscp import solve_rb, total_profile, response, run_field
plt.rcParams.update({'mathtext.fontset': 'dejavuserif', 'font.size': 9, 'legend.frameon': False})
cols = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple']
# ---- Figure 2
fig, axs = plt.subplots(1, 2, figsize=(6.9, 2.8))
for (k, (tau, K, m, DE, DI)), c in zip(SCEN.items(), cols):
    axs[0].plot(simulate(tau, K, m, DE, DI, 151)['NX'], color=c, lw=1.3, label=f'scenario {k}')
axs[0].set_title('(a) New exposures, scenarios of Table 2'); axs[0].set_xlabel('step after seeding'); axs[0].set_ylabel(r'new exposures $p_N(t)$'); axs[0].legend(); axs[0].set_xlim(0, 150)
tau, K, m, DE, DI = SCEN[3]; q = 1 - m; L = DE + DI; T = 151
r = simulate(tau, K, m, DE, DI, T); N = np.concatenate([np.zeros(L + 1), r['NX']])
E = np.array([sum(q ** j * N[t + L + 1 - j] for j in range(1, DE + 1)) for t in range(T)])
I = r['IX']; S = r['S']; Rm = 1 - S - E - I
for y, lab, c, ls in ((S, '$p_S$', '#0072B2', '-'), (E, '$p_E$', '#E69F00', '-'), (I, '$p_I$', '#D55E00', '-'), (Rm, '$p_R$', '#CC79A7', '-'), (10 * r['NX'], '$10\\,p_N$', 'tab:green', '--')):
    axs[1].plot(y, color=c, ls=ls, lw=1.3, label=lab)
axs[1].set_title('(b) Compartments, scenario 3'); axs[1].set_xlabel('step after seeding'); axs[1].set_ylabel('fraction of the unit'); axs[1].legend(ncol=2); axs[1].set_xlim(0, 150)
savefig(fig, 'fig2.png')
# ---- Figure 3
fig, axs = plt.subplots(1, 2, figsize=(6.27, 2.9), sharey=True)
for ax, f, lab in zip(axs, (1.5, 2.0), ('(a) $\\tau_Y = 1.5\\,\\tau_X$', '(b) $\\tau_Y = 2\\,\\tau_X$')):
    for (k, (tau, K, m, DE, DI)), c in zip(SCEN.items(), cols):
        r = simulate(tau, K, m, DE, DI, 241, tauY=f * tau, TsY=0, start='fp')
        ax.plot(r['NY'], color=c, lw=1.3, label=f'scenario {k}'); ax.plot(r['NX'], color=c, lw=1.1, ls=':'); ax.plot(r['NX'] + r['NY'], color=c, lw=1.0, ls='--')
    ax.set_title(lab); ax.set_xlim(0, 240); ax.set_ylim(0, 0.037); ax.set_xticks(range(0, 241, 40)); ax.set_xlabel('step after seeding of Y')
axs[0].set_ylabel('new exposures per step')
fig.legend(*axs[0].get_legend_handles_labels(), loc='upper center', ncol=5, fontsize=8)
axs[1].legend([Line2D([], [], color='0.3', lw=1.3), Line2D([], [], color='0.3', lw=1.1, ls=':'), Line2D([], [], color='0.3', lw=1.0, ls='--')], ['$p_{N,Y}$', '$p_{N,X}$', '$p_{N,X}+p_{N,Y}$'], loc='lower right', fontsize=7.5)
fig.tight_layout(rect=(0, 0, 1, 0.92)); fig.savefig(os.path.join(FIG, 'fig3.png'), dpi=300); plt.close(fig); print('wrote fig3.png')
# ---- Figures 4 and 5
tau, K, m, DE, DI = SCEN[1]; alpha = 0.1
r2 = solve_rb(7, DE, DI, 4.0); r3 = solve_rb(5, DE, DI, 3.0)
prof = {'x1': total_profile(1, 0.0, 0.225, DE, DI, 3000), 'x2': total_profile(7, r2, 0.0375, DE, DI, 3000), 'x3': total_profile(5, r3, 0.7 / 12, DE, DI, 3000)}
names = {'x1': 'CSCP$_{x1}$', 'x2': 'CSCP$_{x2}$', 'x3': 'CSCP$_{x3}$'}
r = simulate(tau, K, m, DE, DI, 601, tauY=2 * tau, TsY=400, start='naive', TsX=0)
fig, ax = plt.subplots(figsize=(5.2, 3.1))
ax.plot(r['NX'] + r['NY'], color='k', lw=1.3, label='transmission unit')
for k, c in zip(prof, cols):
    ax.plot(response(alpha * (r['IX'] + r['IY']), prof[k]['g']), color=c, lw=1.3, label=names[k])
ax.axvline(400, color='0.6', lw=0.8, ls='--'); ax.set_xlim(0, 600); ax.set_ylim(0, None)
ax.set_xlabel('step after seeding of X'); ax.set_ylabel(r'new exposures $p_N$ (unit), $\mathcal{I}$ (CSCPs)'); ax.legend(loc='upper left')
savefig(fig, 'fig4.png')
_, _, IX0 = fixed_point(tau, K, m, DE, DI); pre = 3000
rr = simulate(tau, K, m, DE, DI, 361, tauY=2 * tau, TsY=50, start='fp')   # to step 360: the 1% crossings of X are at 235, 259, 349
fig, ax = plt.subplots(figsize=(5.2, 3.1))
for k, c in zip(prof, cols):
    kx = np.concatenate([np.full(pre, alpha * IX0), alpha * rr['IX']]); ky = np.concatenate([np.zeros(pre), alpha * rr['IY']])
    IXc = response(kx, prof[k]['g'])[pre:]; IYc = response(ky, prof[k]['g'])[pre:]
    ax.plot(IXc + IYc, color=c, lw=1.3, label=names[k]); ax.plot(IYc, color=c, lw=1.3, ls='--'); ax.plot(IXc, color=c, lw=1.3, ls=':')
ax.axvline(50, color='0.6', lw=0.8, ls='--'); ax.set_xlim(0, 360); ax.set_ylim(0, None); ax.set_xlabel('step'); ax.set_ylabel(r'exposure intensity $\mathcal{I}(t)$')
l1 = ax.legend(loc='upper right'); ax.add_artist(l1)
ax.legend([Line2D([], [], color='0.3', lw=1.3), Line2D([], [], color='0.3', lw=1.3, ls='--'), Line2D([], [], color='0.3', lw=1.3, ls=':')], ['X + Y', 'Y', 'X'], loc='center right', bbox_to_anchor=(1, 0.42))
savefig(fig, 'fig5.png')
# ---- Figure 7
KAPPA = 0.025; burn = 3000
fig, ax = plt.subplots(figsize=(4.73, 2.83))
ax.axvspan(200, 220, color='0.85', lw=0, label='perturbation (steps 200–219)'); ax.axhline(0.25, color='0.4', lw=0.8, ls='--')
for (k, n, rb, rho), c in zip((('x1', 1, 0.0, 0.225), ('x2', 7, r2, 0.0375), ('x3', 5, r3, 0.7 / 12)), cols):
    p = run_field(n, DE, DI, rb, rho, lambda t: KAPPA, burn + 600, pert=(burn + 200, burn + 220, 1.5, 1.25))[burn:]
    ax.plot(np.arange(180, 501), p[180:501], color=c, lw=1.3, label=names[k])
h, l = ax.get_legend_handles_labels(); ax.legend(h[1:] + h[:1], l[1:] + l[:1], fontsize=8)
ax.set_xlim(180, 500); ax.set_ylim(0.23, 0.56); ax.set_xlabel('step'); ax.set_ylabel(r'exposure intensity $\mathcal{I}(t)$')
savefig(fig, 'fig7.png')
