import copy
import itertools
import json
import unittest
from fractions import Fraction as F
from src.verification.matrix import (MatrixEnclosure as Box, certify, replay,
    norm_upper, sqrt_upper, inverse_and_det, periodic_return_enclosure)


class ValidatedMatrices(unittest.TestCase):
    def test_exact_diagonal_complete_workflow(self):
        c = certify(Box.point([[2, 0], [0, 3]]), candidate=[1, 0], rhs=[2, 0])
        self.assertEqual(F(c['reference_gap_lower']), 2)
        self.assertEqual(F(c['solution_error_upper']), 0)
        self.assertEqual(F(c['reference_determinant']), 6)
        self.assertTrue(replay(json.loads(json.dumps(c))))

    def test_sharp_scalar(self):
        c = certify(Box([['3/4']], [[1]]), [[1]], candidate=[1], rhs=[1])
        self.assertEqual(F(c['singular_gap_lower']), F(3, 4))
        self.assertEqual(F(c['solution_error_upper']), F(1, 3))

    def test_nonnormal_eigenvalues_are_not_gaps(self):
        # Both eigenvalues are 1; the singular gap is below 1/100.
        c = certify(Box.point([[1, 100], [0, 1]]))
        self.assertLess(F(c['reference_gap_lower']), F(1, 100))
        self.assertGreaterEqual(F(c['reference_gap_lower']), F(1, 101))
        with self.assertRaises(ValueError):
            certify(Box([[1, 100], [0, 1]], [[1, 100], ['1/10', 1]]),
                    [[1, 100], [0, 1]])

    def test_inverse_independent_identity(self):
        for a in ([[0, 2], [3, 1]], [[1, 5, 7], [0, 2, 3], [4, 1, 8]]):
            inv, _ = inverse_and_det(a)
            for i in range(len(a)):
                for j in range(len(a)):
                    self.assertEqual(sum(a[i][k]*inv[k][j] for k in range(len(a))), int(i == j))

    def test_sqrt_and_extreme_rationals(self):
        for q in (F(0), F(2), F(1, 10**100), F(10**1000)):
            u = sqrt_upper(q)
            self.assertGreaterEqual(u*u, q)
            if u:
                self.assertLess((u-F(1, 2**64))**2, q)
        self.assertTrue(replay(certify(Box.point([[F(10**1000)]]))))

    def test_interval_products_all_vertices(self):
        a = Box([[1, -2], [0, 1]], [[2, -1], [1, 2]])
        product = a @ a
        corners = [list(bits) for bits in itertools.product((0, 1), repeat=4)]
        for x in corners:
            for y in corners:
                aa = [[(a.lower, a.upper)[x[2*i+j]][i][j] for j in range(2)] for i in range(2)]
                bb = [[(a.lower, a.upper)[y[2*i+j]][i][j] for j in range(2)] for i in range(2)]
                for i in range(2):
                    for j in range(2):
                        v = sum(aa[i][k]*bb[k][j] for k in range(2))
                        self.assertLessEqual(product.lower[i][j], v)
                        self.assertGreaterEqual(product.upper[i][j], v)

    def test_periodic_order_and_uncertainty(self):
        a = Box.point([[1, 2], [0, 1]])
        b = Box.point([[2, 0], [3, 1]])
        c = periodic_return_enclosure([a, b])
        self.assertEqual(c.lower, ((-1, -4), (-3, -6)))
        uncertain = Box([['49/100', 0], [0, '49/100']],
                        [['51/100', 0], [0, '51/100']])
        receipt = certify(periodic_return_enclosure([uncertain, uncertain]))
        self.assertTrue(replay(receipt))
        self.assertGreater(F(receipt['singular_gap_lower']), 0)

    def test_replay_tampering(self):
        c = certify(Box.point([[1]]), candidate=[1], rhs=[1])
        for key in c:
            changed = copy.deepcopy(c)
            changed[key] = None
            self.assertFalse(replay(changed), key)
        changed = copy.deepcopy(c)
        changed['version'] = True
        self.assertFalse(replay(changed))
        self.assertFalse(replay([]))
        c['unsupported_claim'] = True
        self.assertFalse(replay(c))

    def test_invalid_inputs(self):
        for lo, hi in (([[2]], [[1]]), ([], []), ([[1, 2]], [[1]])):
            with self.assertRaises(ValueError): Box(lo, hi)
        for a in ([[0]], [[1, 2]], [[1.0]], [[True]]):
            with self.assertRaises(ValueError): certify(Box.point(a))
        with self.assertRaises(ValueError): certify(Box.point([[1]]), candidate=[1])
        with self.assertRaises(ValueError): periodic_return_enclosure([])
