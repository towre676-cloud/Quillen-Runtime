"""Quantitative log-determinant variation from derivative and gap bounds."""
import math

def logdet_variation_bound(dimension, derivative_norm_integral, singular_gap):
    """Caller must certify the derivative integral and uniform singular gap."""
    if dimension < 1 or not math.isfinite(derivative_norm_integral) or derivative_norm_integral < 0 or not math.isfinite(singular_gap) or singular_gap <= 0:
        raise ValueError('positive dimension/gap and finite nonnegative integral required')
    return 2 * dimension * derivative_norm_integral / singular_gap
