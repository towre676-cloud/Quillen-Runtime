"""Evaluate spectral and solution-error bounds from supplied norm enclosures."""
import math
from fractions import Fraction


def _nonnegative_finite(value, name):
    if value < 0 or (not isinstance(value, (int, Fraction)) and not math.isfinite(value)):
        raise ValueError(f'{name} must be finite and nonnegative')


def _dimension(dimension):
    if isinstance(dimension, bool) or not isinstance(dimension, int) or dimension < 1:
        raise ValueError('dimension must be a positive integer')


def singular_gap_bound(reference_gap, operator_error=0):
    """Lower singular gap given a reference lower gap and upper operator error.

    Requires operator_error < reference_gap, so the returned margin is positive.
    """
    _nonnegative_finite(reference_gap, 'reference gap')
    _nonnegative_finite(operator_error, 'operator error')
    if operator_error >= reference_gap:
        raise ValueError('operator error must be smaller than the reference gap')
    return reference_gap - operator_error


def logdet_variation_bound(dimension, derivative_norm_integral, singular_gap):
    """Path bound with a uniform lower singular gap and derivative-norm integral."""
    _dimension(dimension)
    _nonnegative_finite(derivative_norm_integral, 'derivative norm integral')
    gap = singular_gap_bound(singular_gap)
    return 2 * dimension * derivative_norm_integral / gap


def logdet_perturbation_bound(dimension, reference_gap, operator_error):
    """Absolute log-det(D*D) error bound, using log1p near zero perturbation."""
    _dimension(dimension)
    singular_gap_bound(reference_gap, operator_error)
    ratio = operator_error / reference_gap
    return -2 * dimension * math.log1p(-ratio)


def logdet_perturbation_rational_bound(dimension, reference_gap, operator_error):
    """Rational upper bound 2*n*epsilon/(s-epsilon); exact for Fraction inputs."""
    _dimension(dimension)
    gap = singular_gap_bound(reference_gap, operator_error)
    return 2 * dimension * operator_error / gap


def solution_error_bound(residual_norm, state_norm, reference_gap, operator_error=0):
    """Bound ||xhat - D^-1 b|| from the reference residual Dhat*xhat - b."""
    _nonnegative_finite(residual_norm, 'residual norm')
    _nonnegative_finite(state_norm, 'state norm')
    gap = singular_gap_bound(reference_gap, operator_error)
    return (residual_norm + operator_error * state_norm) / gap
