"""Rigorous interval arithmetic on Python's decimal module (standard library only).

Every operation is evaluated twice, once in a context that rounds toward minus infinity and
once in a context that rounds toward plus infinity, so the result interval encloses the exact
result of the operation on any pair of points of the operand intervals.  decimal arithmetic is
correctly rounded in the requested direction for +, -, *, / (IEEE 854 / General Decimal
Arithmetic).  Square roots are computed with extra precision and widened by more than the
rounding error.  Decimal input such as '0.035' is represented exactly.
"""
from decimal import Decimal, Context, ROUND_FLOOR, ROUND_CEILING, ROUND_HALF_EVEN

PREC = 80
DN = Context(prec=PREC, rounding=ROUND_FLOOR)
UP = Context(prec=PREC, rounding=ROUND_CEILING)


def set_precision(p):
    global PREC
    PREC = p
    DN.prec = p
    UP.prec = p


class I:
    """Closed interval [lo, hi] with Decimal end points."""
    __slots__ = ('lo', 'hi')

    def __init__(self, lo, hi=None):
        lo = lo if isinstance(lo, Decimal) else Decimal(str(lo))
        hi = lo if hi is None else (hi if isinstance(hi, Decimal) else Decimal(str(hi)))
        if lo > hi:
            raise ValueError('empty interval')
        self.lo, self.hi = lo, hi

    @staticmethod
    def _c(x):
        return x if isinstance(x, I) else I(x)

    def __add__(a, b):
        b = I._c(b); return I(DN.add(a.lo, b.lo), UP.add(a.hi, b.hi))
    __radd__ = __add__

    def __sub__(a, b):
        b = I._c(b); return I(DN.subtract(a.lo, b.hi), UP.subtract(a.hi, b.lo))

    def __rsub__(a, b):
        return I._c(b) - a

    def __neg__(a):
        return I(-a.hi, -a.lo)

    def __mul__(a, b):
        b = I._c(b)
        lo = min(DN.multiply(x, y) for x in (a.lo, a.hi) for y in (b.lo, b.hi))
        hi = max(UP.multiply(x, y) for x in (a.lo, a.hi) for y in (b.lo, b.hi))
        return I(lo, hi)
    __rmul__ = __mul__

    def __truediv__(a, b):
        b = I._c(b)
        if b.lo <= 0 <= b.hi:
            raise ZeroDivisionError('divisor interval contains 0')
        lo = min(DN.divide(x, y) for x in (a.lo, a.hi) for y in (b.lo, b.hi))
        hi = max(UP.divide(x, y) for x in (a.lo, a.hi) for y in (b.lo, b.hi))
        return I(lo, hi)

    def pow_int(a, n):
        """a**n for an integer n >= 0 and a >= 0."""
        if a.lo < 0:
            raise ValueError('pow_int implemented for nonnegative intervals')
        r = I(1); base = a; k = int(n)
        while k:
            if k & 1: r = r * base
            k >>= 1
            if k: base = base * base
        return r

    def sqrt(a):
        if a.lo < 0:
            raise ValueError('sqrt of negative interval')
        ctx = Context(prec=PREC + 10, rounding=ROUND_HALF_EVEN)
        def bounds(x):
            r = ctx.sqrt(x)
            if r == 0: return Decimal(0), Decimal(0)
            eps = Decimal(1).scaleb(r.adjusted() - (PREC + 10) + 2)   # > rounding error
            return DN.subtract(r, eps), UP.add(r, eps)
        lo = bounds(a.lo)[0]; hi = bounds(a.hi)[1]
        return I(max(lo, Decimal(0)), hi)

    def pow_real(a, K):
        """a**K for K an integer or a half-integer (exact Decimal), a > 0."""
        K = K if isinstance(K, Decimal) else Decimal(str(K))
        n = int(K.to_integral_value(rounding=ROUND_FLOOR))
        frac = K - n
        if frac == 0: return a.pow_int(n)
        if frac == Decimal('0.5'): return a.pow_int(n) * a.sqrt()
        raise ValueError('exponent must be an integer or a half-integer')

    def width(a):
        return UP.subtract(a.hi, a.lo)

    def mid(a):
        ctx = Context(prec=PREC + 10, rounding=ROUND_HALF_EVEN)
        return ctx.divide(ctx.add(a.lo, a.hi), 2)

    def mag(a):
        return max(abs(a.lo), abs(a.hi))

    def hull(a, b):
        return I(min(a.lo, b.lo), max(a.hi, b.hi))

    def subset(a, b):
        return b.lo <= a.lo and a.hi <= b.hi

    def __repr__(self):
        return f'[{fmt_dn(self.lo, 20)}, {fmt_up(self.hi, 20)}]'


def fmt_dn(x, digits=25):
    if x == 0: return '0'
    return format(Context(prec=digits, rounding=ROUND_FLOOR).plus(x), 'E')


def fmt_up(x, digits=25):
    if x == 0: return '0'
    return format(Context(prec=digits, rounding=ROUND_CEILING).plus(x), 'E')
