"""Exact real Euclidean matrix-enclosure certificates (standard library only)."""
from fractions import Fraction as F
from math import isqrt
import json
from .check import (singular_gap_bound, logdet_perturbation_rational_bound,
                    solution_error_bound)


def rational(x):
    """No implicit float/decimal rounding; encode measurements as rational strings.

    To enclose a binary float exactly, explicitly pass Fraction.from_float(x).
    """
    if isinstance(x, bool) or not isinstance(x, (int, str, F)):
        raise ValueError('expected integer, rational string, or Fraction')
    return F(x)


def matrix(a):
    a = tuple(tuple(rational(x) for x in row) for row in a)
    if not a or not a[0] or any(len(row) != len(a[0]) for row in a):
        raise ValueError('nonempty rectangular matrix required')
    return a


def sqrt_upper(q):
    """Upper square root on a 2^-64 grid, using integer arithmetic only."""
    q = rational(q)
    if q < 0:
        raise ValueError('negative square')
    scale = 1 << 64
    k = isqrt(q.numerator * scale * scale // q.denominator)
    if F(k*k, scale*scale) < q:
        k += 1
    return F(k, scale)


def norm_upper(a):
    """min(Frobenius, sqrt(||A||_1 ||A||_infinity)), rounded upward."""
    a = matrix(a)
    frob2 = sum(x*x for row in a for x in row)
    rows = max(sum(abs(x) for x in row) for row in a)
    cols = max(sum(abs(x) for x in col) for col in zip(*a))
    return sqrt_upper(min(frob2, rows*cols))


def inverse_and_det(a):
    a = matrix(a)
    n = len(a)
    if len(a[0]) != n:
        raise ValueError('square matrix required')
    work = [list(row)+[F(i == j) for j in range(n)] for i, row in enumerate(a)]
    det = F(1)
    for k in range(n):
        pivot = next((i for i in range(k, n) if work[i][k]), None)
        if pivot is None:
            raise ValueError('singular reference matrix')
        if pivot != k:
            work[k], work[pivot] = work[pivot], work[k]
            det = -det
        p = work[k][k]
        det *= p
        work[k] = [x/p for x in work[k]]
        for i in range(n):
            if i != k:
                p = work[i][k]
                work[i] = [x-p*y for x, y in zip(work[i], work[k])]
    return matrix([row[n:] for row in work]), det


class MatrixEnclosure:
    """Closed entrywise rational intervals; containment is an input hypothesis."""
    def __init__(self, lower, upper):
        self.lower, self.upper = matrix(lower), matrix(upper)
        if (len(self.lower), len(self.lower[0])) != (len(self.upper), len(self.upper[0])):
            raise ValueError('endpoint shape mismatch')
        if any(l > u for lr, ur in zip(self.lower, self.upper) for l, u in zip(lr, ur)):
            raise ValueError('reversed interval')

    @classmethod
    def point(cls, a):
        return cls(a, a)

    def __matmul__(self, other):
        if len(self.lower[0]) != len(other.lower):
            raise ValueError('product shape mismatch')
        lo, hi = [], []
        for i in range(len(self.lower)):
            lr, ur = [], []
            for j in range(len(other.lower[0])):
                lsum = usum = F(0)
                for k in range(len(other.lower)):
                    products = [x*y for x in (self.lower[i][k], self.upper[i][k])
                                for y in (other.lower[k][j], other.upper[k][j])]
                    lsum += min(products)
                    usum += max(products)
                lr.append(lsum)
                ur.append(usum)
            lo.append(lr)
            hi.append(ur)
        return MatrixEnclosure(lo, hi)


def periodic_return_enclosure(steps):
    """I - S_last ... S_first from validated step-map enclosures.

    Includes supplied integration uncertainty via exact interval products.
    Does not itself validate an ODE integrator or an orbit/phase condition.
    """
    steps = list(steps)
    if not steps:
        raise ValueError('at least one step required')
    n = len(steps[0].lower)
    ident = [[int(i == j) for j in range(n)] for i in range(n)]
    product = MatrixEnclosure.point(ident)
    for step in steps:
        if len(step.lower) != n or len(step.lower[0]) != n:
            raise ValueError('square step maps of common dimension required')
        product = step @ product
    return MatrixEnclosure(
        [[ident[i][j]-product.upper[i][j] for j in range(n)] for i in range(n)],
        [[ident[i][j]-product.lower[i][j] for j in range(n)] for i in range(n)])


def certify(enclosure, reference=None, *, candidate=None, rhs=None):
    """Generate JSON-safe, entirely rational certificate for every matrix in a box."""
    lo, hi = enclosure.lower, enclosure.upper
    n = len(lo)
    if len(lo[0]) != n:
        raise ValueError('square enclosure required')
    ref = matrix(reference) if reference is not None else matrix(
        [[(l+u)/2 for l, u in zip(lr, ur)] for lr, ur in zip(lo, hi)])
    if len(ref) != n or len(ref[0]) != n:
        raise ValueError('reference shape mismatch')
    inv, det = inverse_and_det(ref)
    inv_norm = norm_upper(inv)
    s = 1/inv_norm
    radius = [[max(abs(lo[i][j]-ref[i][j]), abs(hi[i][j]-ref[i][j]))
               for j in range(n)] for i in range(n)]
    eps = norm_upper(radius)
    gap = singular_gap_bound(s, eps)
    encode = lambda a: [[str(x) for x in row] for row in a]
    result = dict(version=1, field='real', metric='Euclidean',
                  lower=encode(lo), upper=encode(hi), reference=encode(ref),
                  inverse=encode(inv), reference_determinant=str(det),
                  determinant_sign=1 if det > 0 else -1,
                  inverse_norm_upper=str(inv_norm), reference_gap_lower=str(s),
                  operator_error_upper=str(eps), singular_gap_lower=str(gap),
                  absolute_logdet_error_upper=str(logdet_perturbation_rational_bound(n, s, eps)))
    if (candidate is None) != (rhs is None):
        raise ValueError('candidate and rhs must be supplied together')
    if candidate is not None:
        x, b = tuple(map(rational, candidate)), tuple(map(rational, rhs))
        if len(x) != n or len(b) != n:
            raise ValueError('vector dimension mismatch')
        residual = [[sum(ref[i][j]*x[j] for j in range(n))-b[i]] for i in range(n)]
        eta, state = norm_upper(residual), norm_upper([[v] for v in x])
        result.update(candidate=list(map(str, x)), rhs=list(map(str, b)),
                      residual_norm_upper=str(eta), state_norm_upper=str(state),
                      solution_error_upper=str(solution_error_bound(eta, state, s, eps)))
    return result


def replay(receipt):
    """Recompute all evidence and claims; reject missing, extra or modified fields."""
    if not isinstance(receipt, dict):
        return False
    try:
        computed = certify(MatrixEnclosure(receipt['lower'], receipt['upper']),
                           receipt['reference'], candidate=receipt.get('candidate'),
                           rhs=receipt.get('rhs'))
        return json.dumps(computed, sort_keys=True) == json.dumps(receipt, sort_keys=True)
    except (KeyError, TypeError, ValueError, ZeroDivisionError):
        return False
