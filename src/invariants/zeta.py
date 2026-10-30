"""Finite positive-spectrum log determinant; no zeta/UV extrapolation."""
import math

def positive_logdet(eigenvalues, *, omit_zeros=False):
    values = list(eigenvalues)
    if any(not math.isfinite(x) or x < 0 for x in values):
        raise ValueError('finite nonnegative spectrum required')
    if not omit_zeros and any(x == 0 for x in values):
        raise ValueError('singular: request an explicit primed convention')
    return math.fsum(math.log(x) for x in values if x > 0)
