"""Floating-point simulation of turnover-renewed transmission units (TRTUs) and aggregates (TRTAs).

Notation as in the article: q = 1 - m, d = D_E + 1, L = D_E + D_I, seed 0.01 m at the seeding step.
  p_S(t+1) = q (p_S(t) - p_N(t)) + m
  p_I(t)   = sum_{j=d}^{L} q^j p_N(t-j)
  p_N(t)   = p_S(t) [1 - (1 - tau p_I(t))^K]
Two strains share the escape probability (1 - w)^K, w = tau_X p_I,X + tau_Y p_I,Y, and exposures are
allocated in proportion to tau_V p_I,V (first exposure counts, complete cross-immunity).
"""
import numpy as np


def fixed_point(tau, K, m, DE, DI):
    """Endemic fixed point (p_N*, p_S*, p_I*) by bisection; (0, 1, 0) if there is none."""
    q = 1 - m; d = DE + 1; L = DE + DI
    A = sum(q ** j for j in range(d, L + 1))
    f = lambda N: (1 - q * N / m) * (1 - (1 - tau * A * N) ** K) - N
    lo, hi = 1e-15, m
    if f(lo) <= 0:
        return 0.0, 1.0, 0.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0: lo = mid
        else: hi = mid
    N = 0.5 * (lo + hi)
    return N, 1 - q * N / m, A * N


def simulate(tauX, K, m, DE, DI, T, tauY=None, TsY=None, start='naive', TsX=0, K_of_t=None):
    """Simulate a (one- or two-strain) TRTU for T steps.
    start='naive': fully susceptible, X seeded with 0.01 m at step TsX.
    start='fp'   : constant X-only fixed-point history before step 0.
    Y (if tauY is given) is seeded with 0.01 m at step TsY.  K_of_t(t) optionally overrides K.
    Returns arrays indexed by step t = 0..T-1 (NX, NY, IX, IY, S)."""
    q = 1 - m; d = DE + 1; L = DE + DI; P = L + 5
    NX = np.zeros(T + P); NY = np.zeros(T + P); S = np.ones(T + P)
    IX = np.zeros(T + P); IY = np.zeros(T + P)
    if start == 'fp':
        N0, S0, _ = fixed_point(tauX, K, m, DE, DI)
        NX[:P] = N0; S[:P] = S0
    for i in range(P, T + P):
        t = i - P
        S[i] = (S[i - 1] - NX[i - 1] - NY[i - 1]) * q + m
        IX[i] = sum(q ** j * NX[i - j] for j in range(d, L + 1))
        IY[i] = sum(q ** j * NY[i - j] for j in range(d, L + 1))
        w = tauX * IX[i] + (tauY * IY[i] if tauY else 0.0)
        Kt = K if K_of_t is None else K_of_t(t)
        lam = 1 - (1 - w) ** Kt
        nx = lam * S[i] * tauX * IX[i] / w if w > 0 else 0.0
        ny = lam * S[i] * tauY * IY[i] / w if (w > 0 and tauY) else 0.0
        sx = 0.01 * m if (start == 'naive' and t == TsX) else 0.0
        sy = 0.01 * m if (TsY is not None and t == TsY) else 0.0
        if sx + sy > 0:
            nx = nx * (1 - (sx + sy) / S[i]) + sx
            ny = ny * (1 - (sx + sy) / S[i]) + sy
        NX[i] = nx; NY[i] = ny
    sl = slice(P, T + P)
    return dict(NX=NX[sl], NY=NY[sl], IX=IX[sl], IY=IY[sl], S=S[sl])


def ramp(tau, K0, K1, n, m, DE, DI, T=150):
    """Start at the fixed point for K0; at step s >= 1 the contact intensity is K0 + (K1-K0) min(s,n)/n.
    Returns p_N(s) for s = 1..T (index s-1)."""
    q = 1 - m; d = DE + 1; L = DE + DI; P = L + 5
    N0, S0, _ = fixed_point(tau, K0, m, DE, DI)
    N = np.full(T + P, N0); S = np.full(T + P, S0)
    for i in range(P, T + P):
        s = i - P + 1
        K = K0 + (K1 - K0) * min(s, n) / n
        S[i] = (S[i - 1] - N[i - 1]) * q + m
        I = sum(q ** j * N[i - j] for j in range(d, L + 1))
        N[i] = S[i] * (1 - (1 - tau * I) ** K)
    return N[P:]


def aggregate(x, weights, offsets, T, pre=None):
    """Weighted sum of shifted copies of a unit trajectory x (index = steps after the unit's seeding)."""
    a = np.zeros(T)
    for w, o in zip(weights, offsets):
        fill = x[0] if pre is None else pre
        a += w * np.concatenate([np.full(o, fill), x[:T - o]])
    return a


def successive_extrema(x, kinds, start=1, horizon=None):
    """Successive global extrema after each other, e.g. ('max','min','max'); returns [(value, index)]."""
    out = []; s = start; e = len(x) if horizon is None else horizon
    for k in kinds:
        seg = x[s:e]
        i = s + int(np.argmax(seg) if k == 'max' else np.argmin(seg))
        out.append((float(x[i]), i)); s = i + 1
    return out


def replacement_times(NX, NY, X_star, TsY):
    """Table S3 / Table 2 conventions.
    counted_1_95: steps from the first step after seeding at which the Y share reaches 1% to the first
                  step at which it reaches 95%, both included;
    elapsed_to_1: elapsed steps after the seeding of Y until X <= 1% of its X-only level;
    elapsed_90_1: elapsed steps from the first step with X <= 90% to the first step with X <= 1%."""
    t = np.arange(len(NX))
    share = NY / np.where(NX + NY > 0, NX + NY, 1)
    t1 = int(np.argmax((share >= 0.01) & (t > TsY)))
    t95 = int(np.argmax((share >= 0.95) & (t > TsY)))
    tx1 = int(np.argmax((NX <= 0.01 * X_star) & (t > TsY)))
    tx90 = int(np.argmax((NX <= 0.90 * X_star) & (t > TsY)))
    return dict(counted_1_95=t95 - t1 + 1, elapsed_to_1=tx1 - TsY, elapsed_90_1=tx1 - tx90)
