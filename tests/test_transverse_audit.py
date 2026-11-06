"""Additional integration regressions from the November 6 independent audit."""
import copy
from fractions import Fraction as F
import importlib.util
import unittest

from examples.check_transverse import check
from src.verification.matrix import MatrixEnclosure as Box
from src.verification.transverse import exact_return_certificate, replay, sampled_return_certificate
from src.verification.dynacomp_bridge import sampled_step

Q = [[1], [0]]
L = [[1, 0]]


class ExactAudit(unittest.TestCase):
    def test_product_invariance_does_not_override_segment_checks(self):
        # The total map is (1/2) I, but neither factor preserves span(e1).
        a = [[0, 1], [1, 0]]
        b = [[0, '1/2'], ['1/2', 0]]
        with self.assertRaises(ValueError):
            exact_return_certificate([a, b], Q, L)
        self.assertTrue(check(exact_return_certificate([[['1/2', 0], [0, '1/2']]], Q, L)))

    def test_neutral_jordan_chain_cannot_be_declared_fixed(self):
        s = [['1/2', 0, 0], [0, 1, 1], [0, 0, 1]]
        q, l = [[1], [0], [0]], [[1, 0, 0]]
        c = exact_return_certificate([s], q, l)
        self.assertEqual(c['reduction'], 'invariant-subspace-only')
        self.assertTrue(check(c))
        with self.assertRaises(ValueError):
            exact_return_certificate([s], q, l, neutral=[[0, 0], [1, 0], [0, 1]])

    def test_final_only_fixed_neutral_complement_is_unsupported(self):
        # Neutral on the final product, but not pointwise fixed by each step.
        s = [['1/2', 0], [0, -1]]
        self.assertTrue(check(exact_return_certificate([s, s], Q, L)))
        with self.assertRaises(ValueError):
            exact_return_certificate([s, s], Q, L, neutral=[[0], [1]])

    def test_basis_change_preserves_determinant_but_not_coordinate_gap(self):
        s = [['1/2', 0, 0], [0, '1/4', 0], [0, 0, 1]]
        a = exact_return_certificate([s], [[1, 0], [0, 1], [0, 0]],
                                     [[1, 0, 0], [0, 1, 0]], neutral=[[0], [0], [1]])
        b = exact_return_certificate([s], [[1, 10], [0, 1], [0, 0]],
                                     [[1, -10, 0], [0, 1, 0]], neutral=[[0], [0], [1]])
        self.assertEqual(a['certificate']['reference_determinant'],
                         b['certificate']['reference_determinant'])
        self.assertLess(F(b['certificate']['singular_gap_lower']),
                        F(a['certificate']['singular_gap_lower']))
        self.assertTrue(replay(b))
        self.assertTrue(check(b))

    def test_nonnormal_eigenvalue_distance_is_not_claimed_as_gap(self):
        c = exact_return_certificate([[['1/2', 100, 0], [0, '1/2', 0], [0, 0, 1]]],
                                     [[1, 0], [0, 1], [0, 0]], [[1, 0, 0], [0, 1, 0]])
        # Both eigenvalue distances to 1 equal 1/2; the singular gap is tiny.
        self.assertLess(F(c['certificate']['singular_gap_lower']), F(1, 100))
        self.assertTrue(check(c))


@unittest.skipUnless(importlib.util.find_spec('dynacomp'), 'pinned DynaComp required')
class SampledAudit(unittest.TestCase):
    def test_nested_source_and_assumption_tampering_rejected_by_both(self):
        _, step = sampled_step([[-1, 0], [0, 0]], 1)
        c = sampled_return_certificate([step], Q, L, neutral=[[0], [1]])
        for key in step['source_sha256']:
            d = copy.deepcopy(c)
            d['steps'][0]['source_sha256'][key] = '0' * 64
            self.assertFalse(replay(d))
            self.assertFalse(check(d))
        for field in ('nonlinear_orbit_verified', 'external_invariance_assumed'):
            d = copy.deepcopy(c)
            d['steps'][0][field] = True
            self.assertFalse(replay(d))
            self.assertFalse(check(d))

    def test_noncommuting_exponential_against_independent_scalar_series(self):
        a = [[-1, 1, 0], [0, -1, 0], [0, 0, 0]]
        b = [[-1, 0, 0], [1, -1, 0], [0, 0, 0]]
        c = sampled_return_certificate([sampled_step(x, 1)[1] for x in (a, b)],
            [[1, 0], [0, 1], [0, 0]], [[1, 0, 0], [0, 1, 0]], neutral=[[0], [0], [1]])
        term = total = F(1)
        for k in range(1, 41):
            term *= F(-1, k)
            total += term
        exp_lo, exp_hi = total + term * F(-1, 41), total
        e2_lo, e2_hi = exp_lo**2, exp_hi**2
        # exp(B2) exp(B1) = exp(-2) [[1,1],[1,2]].
        for i, row in enumerate([[1, 1], [1, 2]]):
            for j, v in enumerate(row):
                self.assertLessEqual(F(c['certificate']['lower'][i][j]), int(i == j)-v*e2_hi)
                self.assertGreaterEqual(F(c['certificate']['upper'][i][j]), int(i == j)-v*e2_lo)
        self.assertGreater(F(c['certificate']['singular_gap_lower']), F(3, 5))

    def test_singular_mass_and_sampling_limits_fail(self):
        for h, kw in [(0, {}), (-1, {}), (F(9), {}),
                      (1, {'mass': [[1, 0], [0, 0]]})]:
            with self.subTest(h=h, kw=kw), self.assertRaises(ValueError):
                sampled_step([[-1, 0], [0, 0]], h, **kw)
        with self.assertRaises(ValueError):
            sampled_step([[int(i == j) for j in range(12)] for i in range(12)], 1)


if __name__ == '__main__':
    unittest.main()
