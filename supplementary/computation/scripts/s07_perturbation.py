"""Corollary 3.3.6, numerical illustration and Fig. 6: in-cluster probability x1.5 and cross-cluster rate x1.25
during steps 200-219 after a long burn-in at constant seeding kappa = 0.025 (baseline exposure intensity 0.25,
the level behind the values printed in the article).  The field is linear in kappa, so the relative columns
and the step counts do not depend on kappa."""
import numpy as np
from _common import *
from trtu_model.cscp import solve_rb, run_field, cluster_profile
DE, DI = 3, 4
r2 = solve_rb(7, DE, DI, 4.0); r3 = solve_rb(5, DE, DI, 3.0)
systems = [('x1', 1, 0.0, 0.225), ('x2', 7, r2, 0.0375), ('x3', 5, r3, 0.7 / 12)]
KAPPA = 0.025; burn = 3000; rows = []; fig, ax = plt.subplots(figsize=(4.5, 3.0))
for name, n, rb, rho in systems:
    T = burn + 600
    base = run_field(n, DE, DI, rb, rho, lambda t: KAPPA, T)
    per = run_field(n, DE, DI, rb, rho, lambda t: KAPPA, T, pert=(burn + 200, burn + 220, 1.5, 1.25))
    b = base[burn:]; p = per[burn:]; ex = p - b; rel = np.abs(ex) / b
    within = next(t for t in range(220, len(rel)) if np.all(rel[t:] <= 0.01)) - 220
    GCp = cluster_profile(n, 1.5 * rb, DE, DI)[0][-1]
    rows.append([name, round(b[100], 4), round(GCp, 3), round(GCp * 1.25 * rho * DI, 3), round(p.max(), 4), round(p.max() / b[100], 4), int(np.argmax(p)), round(ex[250], 4), round(ex[250] / b[100], 4), within])
    ax.plot(np.arange(150, 500), p[150:500], lw=1.1, label=name)
write_csv('cor336_perturbation.csv', ['system', 'baseline', 'perturbed_in-cluster_size', 'perturbed_cluster_reproduction_number', 'peak', 'peak/baseline', 'step_of_peak', 'excess_at_step_250', 'excess_at_step_250/baseline', 'post-intervention_steps_until_within_1%'], rows)
ax.axvspan(200, 220, color='0.85', lw=0); ax.set_xlabel('step'); ax.set_ylabel('exposure intensity'); ax.legend(frameon=False)
savefig(fig, 'fig6_perturbation.png')
