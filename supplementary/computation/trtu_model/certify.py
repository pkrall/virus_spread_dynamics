"""Computer-assisted proofs for the TRTU examples (Theorem 3.1.3, Table 1, Tables S1 and S2).

Delayed state X_t = (p_S(t), p_N(t-1), ..., p_N(t-L)),  L = D_E + D_I,  d = D_E + 1,  q = 1 - m.
One-step map  F(s, h_1..h_L) = (q s u^K + m,  s (1 - u^K),  h_1, ..., h_{L-1}),
with u = 1 - tau * sum_{j=d}^{L} q^j h_j.  All arithmetic is outward-rounded (interval.py).

Certificate for one trajectory:
  1. fixed point: interval evaluation of H(N) = [1 - (q/m) N](1 - [1 - tau q^d S N]^K) - N at
     two decimal end points shows a sign change; by uniqueness (Lemma 3.1.2) the endemic fixed
     point lies in between.  S = sum_{i<D_I} q^i.
  2. enclosures of the trajectory up to a step T at which the delayed state lies in the box
     B = {x : |x - c|_inf <= r}, c the centre of the fixed-point enclosure.
  3. an enclosure sequence (tube) B_0 = B, B_{k+1} ⊇ F(B_k), k < P, interval Jacobians on the
     B_k and an enclosure of DF^P on B; if its maximum row-sum norm is < 1 (and the fixed-point
     uncertainty is small enough), F^P maps B into itself and contracts there, so the trajectory
     converges and stays in the tube for all t >= T.
  4. extremal values: global maximum / minimum of the enclosed trajectory in the window
     [1, T] with strict separation from all other enclosures and from the tube range.
"""
from decimal import Decimal as D
from .interval import I, DN, UP, fmt_dn, fmt_up


class Unit:
    def __init__(self, tau, K, m, DE, DI):
        self.tau, self.K, self.m = D(str(tau)), D(str(K)), D(str(m))
        self.DE, self.DI = DE, DI
        self.d, self.L = DE + 1, DE + DI
        self.q = I(1) - I(self.m)
        self.qp = [self.q.pow_int(j) for j in range(self.L + 1)]       # q^j as intervals

    # ---------------------------------------------------------------- map and Jacobian
    def u(self, x, tau=None):
        tau = I(self.tau if tau is None else tau)
        s = I(0)
        for j in range(self.d, self.L + 1):
            s = s + self.qp[j] * x[j]
        return I(1) - tau * s

    def F(self, x, K=None):
        K = self.K if K is None else K
        u = self.u(x)
        if u.lo <= 0: raise ValueError('state outside the domain U')
        uk = u.pow_real(K)
        s = x[0]
        return [self.q * s * uk + I(self.m), s * (I(1) - uk)] + x[1:self.L]

    def DF(self, x, K=None):
        K = self.K if K is None else K
        n = self.L + 1
        u = self.u(x); uk = u.pow_real(K); ukm1 = u.pow_real(K - 1)
        s = x[0]; tau = I(self.tau)
        J = [[I(0) for _ in range(n)] for _ in range(n)]
        J[0][0] = self.q * uk
        J[1][0] = I(1) - uk
        for j in range(self.d, self.L + 1):
            g = s * I(K) * ukm1 * tau * self.qp[j]          # derivative of s(1 - u^K) w.r.t. h_j
            J[0][j] = -(self.q * g)
            J[1][j] = g
        for k in range(1, self.L):
            J[k + 1][k] = I(1)
        return J

    # ---------------------------------------------------------------- fixed point
    def H(self, N):
        N = I(N)
        S = I(0)
        for i in range(self.DI):
            S = S + self.qp[i]
        a = I(1) - (self.q / I(self.m)) * N
        b = I(1) - (I(1) - I(self.tau) * self.qp[self.d] * S * N).pow_real(self.K)
        return a * b - N

    def fixed_point(self, N_lo, N_hi):
        h_lo, h_hi = self.H(N_lo), self.H(N_hi)
        ok = h_lo.lo > 0 and h_hi.hi < 0
        N = I(D(N_lo), D(N_hi))
        S = I(1) - (self.q / I(self.m)) * N
        return dict(ok=ok, H_lo=h_lo, H_hi=h_hi, N=N, X=[S] + [N] * self.L)


def refine_root(unit, guess, digits=40):
    """Floating guess -> decimal bracket of width 10^-digits by bisection on point evaluations.
    (Only the bracket is used afterwards; its validity is checked rigorously in fixed_point.)"""
    from decimal import localcontext
    with localcontext() as ctx:
        ctx.prec = 120
        lo, hi = D(str(guess)) * D('0.99'), D(str(guess)) * D('1.01')
        mid_val = lambda x: unit.H(x).mid()
        assert mid_val(lo) > 0 > mid_val(hi), 'bad initial bracket'
        step = D(10) ** (-digits)
        while hi - lo > step:
            mid = (lo + hi) / 2
            if mid_val(mid) > 0: lo = mid
            else: hi = mid
        lo = lo.quantize(step, rounding='ROUND_FLOOR'); hi = hi.quantize(step, rounding='ROUND_CEILING')
    return lo, hi


# -------------------------------------------------------------------- trajectories
def naive_trajectory(unit, T):
    """Seed 0.01 m at step 0 in an epidemiologically naive unit; returns X_1..X_{T+1} and N(0..T)."""
    sig = I(D('0.01')) * I(unit.m)
    X = [unit.q * (I(1) - sig) + I(unit.m), sig] + [I(0)] * (unit.L - 1)
    N = [sig]; Xs = [None, X]
    for t in range(1, T + 1):
        X = unit.F(X); N.append(X[1]); Xs.append(X)
    return Xs, N


def ramp_trajectory(unit, K0, K1, n, X_start, T):
    """Start at the fixed-point state for K0 (X_1 = X*), K(s) = K0 + (K1-K0) min(s,n)/n; N(s), s = 1..T."""
    K0, K1 = D(str(K0)), D(str(K1))
    X = list(X_start); Xs = [None, X]; N = [None]
    for s in range(1, T + 1):
        K = K0 + (K1 - K0) * D(min(s, n)) / D(n)        # exact: K is an integer or a half-integer here
        X = unit.F(X, K); N.append(X[1]); Xs.append(X)
    return Xs, N


def dist_inf(x, c):
    return max(max(UP.subtract(xi.hi, ci), UP.subtract(ci, xi.lo)) for xi, ci in zip(x, c))


def contraction(unit, X_fp, r, P):
    """Tube B_0..B_P around the fixed-point enclosure and an enclosure of DF^P on B_0."""
    c = [xi.mid() for xi in X_fp]
    w = dist_inf(X_fp, c)                      # rigorous bound on |c - X*|_inf
    B = [I(DN.subtract(ci, r), UP.add(ci, r)) for ci in c]
    tube = [B]; M = None
    for k in range(P):
        J = unit.DF(tube[-1])
        M = J if M is None else matmul(J, M)
        tube.append(unit.F(tube[-1]))
    L_bound = max(sum_up([Mij.mag() for Mij in row]) for row in M)
    invariant = UP.add(UP.multiply(L_bound, UP.add(r, w)), w) <= r
    return dict(center=c, B=B, tube=tube, L=L_bound, ok=(L_bound < 1 and invariant), w=w)


def sum_up(xs):
    s = D(0)
    for x in xs: s = UP.add(s, x)
    return s


def matmul(A, B):
    n = len(A)
    out = []
    for i in range(n):
        row = []
        for j in range(n):
            acc = I(0)
            for k in range(n):
                if (A[i][k].lo == 0 and A[i][k].hi == 0) or (B[k][j].lo == 0 and B[k][j].hi == 0):
                    continue
                acc = acc + A[i][k] * B[k][j]
            row.append(acc)
        out.append(row)
    return out


def tail_range(tube):
    """Range of p_N(t) for all t >= T when X_T lies in B: union of the second coordinate of B_1..B_P."""
    r = tube[1][1]
    for b in tube[2:]: r = r.hull(b[1])
    return r


def certify_extrema(N, steps, kinds, tail):
    """kinds e.g. ('max','min','max'): successive global extrema of N over the window after the previous one.
    Each must be strictly separated from every other enclosure in its search range and from the tail range."""
    out = []; start = steps[0]
    for kind in kinds:
        rng = [t for t in steps if t >= start]
        key = (lambda t: N[t].mid()) if kind == 'max' else (lambda t: -N[t].mid())
        t0 = max(rng, key=key)
        others = [N[t] for t in rng if t != t0] + [tail]
        if kind == 'max':
            ok = all(N[t0].lo > o.hi for o in others)
        else:
            ok = all(N[t0].hi < o.lo for o in others)
        out.append(dict(kind=kind, step=t0, enclosure=N[t0], ok=ok))
        start = t0 + 1
    return out
