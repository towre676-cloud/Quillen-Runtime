"""Evaluate the exact Hopf model and its return-map error certificate."""
import json
import math
import sys
from pathlib import Path
from fractions import Fraction as F

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.invariants.zeta import positive_logdet
from src.verification.check import (logdet_perturbation_bound, singular_gap_bound,
                                    solution_error_bound, logdet_perturbation_rational_bound)

rows = []
for r in [.5, 1., 2.]:
    lam = (1-r)**2
    rows.append({
        'r': r,
        'mu': math.log(r)/(2*math.pi),
        'laplacian_eigenvalues': [lam, lam],
        'real_nullity': 2 if r == 1 else 0,
        'logdet': positive_logdet([lam, lam]) if lam else None,
        'primed_logdet': positive_logdet([lam, lam], omit_zeros=True),
    })

certificate = {
    'reference_singular_gap': .5,
    'operator_error_upper': .01,
    'reference_residual_norm_upper': 1e-6,
    'candidate_state_norm_upper': 1.,
    'singular_gap_lower': singular_gap_bound(.5, .01),
    'logarithmic_error_bound_evaluation': logdet_perturbation_bound(2, .5, .01),
    'absolute_logdet_error_upper_exact': str(logdet_perturbation_rational_bound(2, F(1, 2), F(1, 100))),
    'singular_gap_lower_exact': str(singular_gap_bound(F(1, 2), F(1, 100))),
    'solution_error_upper_exact': str(solution_error_bound(F(1, 1000000), F(1), F(1, 2), F(1, 100))),
    'solution_error_upper': solution_error_bound(1e-6, 1., .5, .01),
}
print(json.dumps({
    'model': 'hopf_equilibrium_monodromy',
    'arithmetic': 'binary64',
    'initial_value_spectrum': [(3-math.sqrt(5))/2, (3+math.sqrt(5))/2],
    'periodic_rows': rows,
    'return_map_certificate': certificate,
}, indent=2, allow_nan=False))
