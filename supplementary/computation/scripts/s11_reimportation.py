"""Remark (re-importation of the incumbent) after Lemma A.1.2: a numerical experiment outside the model.

Two-strain transmission unit of Section 2.1 (scenario 1 of Table 2, tau_Y = 2 tau_X), started from the X-only
endemic fixed-point history, Y seeded with 0.01 m at step 0.  Modification: X is seeded at every step
t >= 1 with sigma_X(t) = eps * m (resident seeding rule of Section 2.1, not the admission rule of Section 3.2).
Output: limit of the endogenous X exposures relative to the X-only fixed point, and the first step at which
they fall to 1% of that level."""
from _common import *
from trtu_model.trtu import fixed_point
tau_x, K, m, DE, DI = SCEN[1]; tau_y = 2 * tau_x
q, d, L = 1 - m, DE + 1, DE + DI


def run(eps, T=6000):
    n0, s0, _ = fixed_point(tau_x, K, m, DE, DI)
    off = L + 2
    ntx = [n0] * off; nty = [0.0] * off; s = [s0] * off; nx_end = []
    for step in range(T):
        t = len(s)
        s_t = (s[t - 1] - ntx[t - 1] - nty[t - 1]) * q + m
        ix = sum(q ** j * ntx[t - j] for j in range(d, L + 1))
        iy = sum(q ** j * nty[t - j] for j in range(d, L + 1))
        w = ix * tau_x + iy * tau_y
        lam = 1 - (1 - w) ** K
        nx = lam * s_t * ix * tau_x / w if ix > 0 else 0.0
        ny = lam * s_t * iy * tau_y / w if iy > 0 else 0.0
        sx = eps * m if step >= 1 else 0.0
        sy = 0.01 * m if step == 0 else 0.0
        s.append(s_t); ntx.append(nx * (1 - (sx + sy) / s_t) + sx); nty.append(ny * (1 - (sx + sy) / s_t) + sy)
        nx_end.append(nx)
    first = next((i for i in range(1, T) if nx_end[i] <= 0.01 * n0), None)
    return nx_end[-1] / n0, first


rows = []
for eps in (0.0, 1e-5, 1e-4, 1e-3):
    floor, first = run(eps)
    rows.append([eps, round(100 * floor, 4), first])
write_csv('remark_reimportation.csv', ['eps', 'limit_of_p_N,X_in_%_of_X-only_fixed_point', 'first_step_p_N,X<=1%'], rows)
