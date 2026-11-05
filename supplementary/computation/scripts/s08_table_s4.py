"""Table S4 (Corollary 3.3.8): total seeding k, Y share rising linearly from 0% at step 0 to 100% at step 50.
Right group: exact cascade profiles.  Left group (n = 10, D_I = 20, perfect isolation): the within-cluster
state space is too large for the exact chain, so the age profile is estimated by Monte Carlo (first-passage
percolation with geometric transmission delays, 2,000,000 clusters per strain, seed 1); final sizes are exact
(Reed-Frost)."""
import numpy as np
from _common import *
from trtu_model.cscp import total_profile, reed_frost_mean
rng = np.random.default_rng(1)

def mc_profile(n, rb, DE, DI, N=2_000_000, chunk=200_000, A=400):
    cnt = np.zeros(A + 1)
    for _ in range(N // chunk):
        G = rng.geometric(rb, size=(chunk, n, n)).astype(float); G[G > DI] = np.inf
        W = DE + G; Tm = np.full((chunk, n), np.inf); Tm[:, 0] = 0
        for _ in range(n): Tm = np.minimum(Tm, np.min(Tm[:, :, None] + W, axis=1))
        Tf = Tm[np.isfinite(Tm)].astype(int); cnt += np.bincount(Tf, minlength=A + 1)[:A + 1]
    return cnt / N

def table(gX, gY, label, steps=(0, 10, 25, 40, 50, 60, 80)):
    """Seeding share of Y rises linearly from 0 at step 0 to 1 at step 50; X-only seeding before step 0."""
    A = len(gX); t = np.arange(-A, 200); kY = np.clip(t / 50, 0, 1); kX = 1 - kY
    IX = np.convolve(kX, gX)[:len(t)]; IY = np.convolve(kY, gY)[:len(t)]
    tot = IX + IY; out = []
    for s in steps:
        j = int(np.where(t == s)[0][0]); out.append([label, s, round(100 * kY[j]), round(IX[j], 4), round(IY[j], 4), round(tot[j], 4)])
    m0 = np.where(t >= 0)[0]; jm = m0[np.argmax(tot[m0])]
    out.append([label, f'max (step {t[jm]})', round(100 * kY[jm]), round(IX[jm], 4), round(IY[jm], 4), round(tot[jm], 4)])
    out.append([label, 'stationary', 100, 0.0, round(gY.sum(), 4), round(gY.sum(), 4)])
    return out

PX = total_profile(10, 0.05, 0.002, 3, 5, A=600); PY = total_profile(10, 0.15, 0.0025, 3, 5, A=600)
rows = table(PX['g'][:600], PY['g'][:600], 'right (exact)')
gX = mc_profile(10, 0.05, 3, 20); gY = mc_profile(10, 0.15, 3, 20)
rows += table(gX[:400], gY[:400], 'left (Monte Carlo)')
write_csv('table_S4.csv', ['group', 'step', 'Y_share_of_seeding_%', 'I_X/k', 'I_Y/k', '(I_X+I_Y)/k'], rows)
fin = [['left', 'X', reed_frost_mean(10, 1 - 0.95 ** 20), round(gX.sum(), 5)], ['left', 'Y', reed_frost_mean(10, 1 - 0.85 ** 20), round(gY.sum(), 5)],
       ['right', 'X (within-cluster, total)', round(PX['GammaC'], 8), round(PX['Gamma_tot'], 8)], ['right', 'Y (within-cluster, total)', round(PY['GammaC'], 8), round(PY['Gamma_tot'], 8)],
       ['right', 'cluster reproduction factors X, Y', round(PX['Rc'], 8), round(PY['Rc'], 8)]]
write_csv('table_S4_cascade_sizes.csv', ['group', 'quantity', 'exact', 'Monte_Carlo_or_total'], fin)
