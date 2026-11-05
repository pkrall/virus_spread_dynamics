"""Continuously seeded complementary populations (CSCPs): exact expected cascade profiles.

Within a cluster of size n (complete graph), every infectious member exposes every susceptible member
with probability rb per step (chain binomial); a seeded cluster starts with one exposed member.
D_E exposed steps and D_I infectious steps.  Each infectious individual seeds new clusters at the
expected rate rho per step (cross-cluster transmission).  Expectations are exact: the full
distribution of the within-cluster Markov chain is propagated.
"""
import numpy as np
import scipy.sparse as sp
from math import comb


def cluster_profile(n, rb, DE, DI, horizon=None):
    """Within-cluster expected cumulative cases Gamma_C(a), infectious members i(a) and new exposures g_C(a), a = 0..H."""
    H = horizon if horizon is not None else n * (DE + DI) + 2
    start = (n - 1, 1) + (0,) * (DE - 1) + (0,) * DI + (0,)
    dist = {start: 1.0}
    G = np.zeros(H + 1); ib = np.zeros(H + 1); g = np.zeros(H + 1); q = 1.0 - rb
    for a in range(H + 1):
        for st, p in dist.items():
            G[a] += p * (n - st[0]); ib[a] += p * sum(st[1 + DE:1 + DE + DI]); g[a] += p * st[1]
        if a == H: break
        nd = {}
        for st, p in dist.items():
            S = st[0]; E = st[1:1 + DE]; Iv = st[1 + DE:1 + DE + DI]; R = st[-1]
            pe = 1.0 - q ** sum(Iv)
            for x in range(S + 1):
                if pe == 0.0 and x > 0: break
                pr = comb(S, x) * pe ** x * (1 - pe) ** (S - x)
                if pr == 0.0: continue
                ns = (S - x, x) + E[:-1] + (E[-1],) + Iv[:-1] + (R + Iv[-1],)
                nd[ns] = nd.get(ns, 0.0) + p * pr
        dist = nd
    return G, ib, g


def total_profile(n, rb, rho, DE, DI, A=400):
    """Expected new exposures g(a) by cascade age a for one seeded cluster including cross-cluster seeding."""
    GC, ib, gC = cluster_profile(n, rb, DE, DI)
    HC = len(GC) - 1
    beta = np.zeros(A + 1)
    for a in range(1, A + 1):
        beta[a] = rho * (ib[a - 1] if a - 1 <= HC else 0.0)
    B = np.zeros(A + 1); B[0] = 1.0                 # expected number of clusters newly seeded at cascade age a
    for a in range(1, A + 1):
        B[a] = np.dot(B[:a], beta[a:0:-1])
    gCpad = np.zeros(A + 1); gCpad[:min(A, HC) + 1] = gC[:min(A, HC) + 1]
    g = np.convolve(B, gCpad)[:A + 1]
    GammaC = GC[-1]; Rc = GammaC * rho * DI
    return dict(g=g, Gam=np.cumsum(g), Gamma_tot=GammaC / (1 - Rc) if Rc < 1 else np.inf, GammaC=GammaC, Rc=Rc, B=B)


def solve_rb(n, DE, DI, target, lo=0.0, hi=1.0, iters=80):
    """In-cluster transmission probability giving the expected within-cluster cascade size `target` (bisection)."""
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if cluster_profile(n, mid, DE, DI)[0][-1] < target: lo = mid
        else: hi = mid
    return 0.5 * (lo + hi)


def response(kappa, g):
    """Exposure intensity I(t) = sum_a g(a) kappa(t - a) (Corollary 3.3.4)."""
    return np.convolve(kappa, g)[:len(kappa)]


def reed_frost_mean(n, pi):
    """Expected final size of a Reed-Frost epidemic in a cluster of n started by one case."""
    q = 1 - pi; p = [0, 1.0]
    for k in range(2, n + 1):
        p.append(1 - sum(comb(k - 1, j - 1) * p[j] * q ** (j * (k - j)) for j in range(1, k)))
    return sum(k * comb(n - 1, k - 1) * p[k] * q ** (k * (n - k)) for k in range(1, n + 1))


# ---------------------------------------------------------------- perturbations (Corollary 3.3.6)
def _chain(n, DE, DI, rbs):
    start = (n - 1, 1) + (0,) * (DE - 1) + (0,) * DI + (0,)
    def succ(st, rb):
        S = st[0]; E = st[1:1 + DE]; Iv = st[1 + DE:1 + DE + DI]; R = st[-1]
        pe = 1 - (1 - rb) ** sum(Iv); out = []
        for x in range(S + 1):
            pr = comb(S, x) * pe ** x * (1 - pe) ** (S - x)
            if pr: out.append(((S - x, x) + E[:-1] + (E[-1],) + Iv[:-1] + (R + Iv[-1],), pr))
        return out
    active = lambda st: sum(st[1:-1]) > 0
    states = [start]; idx = {start: 0}; k = 0
    while k < len(states):
        for ns, _ in succ(states[k], 0.5):
            if active(ns) and ns not in idx:
                idx[ns] = len(states); states.append(ns)
        k += 1
    Ps = []
    for rb in rbs:
        rows, cols, vals = [], [], []
        for i, st in enumerate(states):
            for ns, pr in succ(st, rb):
                if ns in idx: rows.append(i); cols.append(idx[ns]); vals.append(pr)
        Ps.append(sp.csr_matrix((vals, (rows, cols)), shape=(len(states), len(states))))
    iC = np.array([sum(st[1 + DE:1 + DE + DI]) for st in states], float)
    eC = np.array([st[1] for st in states], float)
    return Ps, iC, eC


def run_field(n, DE, DI, rb, rho, kappa, T, pert=None):
    """Expected exposures e(t), t = 1..T, of a CSCP under seeding kappa(t) (expected active-cluster field).
    pert = (t0, t1, rb_multiplier, rho_multiplier) applies the multipliers for t0 <= t < t1."""
    rbp = rb * (pert[2] if pert else 1)
    Ps, iC, eC = _chain(n, DE, DI, [rb, rbp])
    phi = np.zeros(Ps[0].shape[0]); out = np.zeros(T + 1)
    for t in range(1, T + 1):
        inP = bool(pert) and pert[0] <= t < pert[1]
        P = Ps[1] if inP else Ps[0]
        r = rho * (pert[3] if inP else 1)
        new = phi @ P
        new[0] += kappa(t) + r * (iC @ phi)
        phi = new; out[t] = eC @ phi
    return out
