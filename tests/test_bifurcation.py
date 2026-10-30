import math
import unittest
from fractions import Fraction as F
from src.core.mapping_cone import defect_matrix
from src.invariants.zeta import positive_logdet
from src.verification.check import logdet_variation_bound


def det(a):
    a = [[F(x) for x in row] for row in a]
    result = F(1)
    for k in range(len(a)):
        pivot = next((i for i in range(k, len(a)) if a[i][k]), None)
        if pivot is None:
            return F(0)
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            result = -result
        p = a[k][k]
        result *= p
        for i in range(k+1, len(a)):
            ratio = a[i][k]/p
            for j in range(k+1, len(a)):
                a[i][j] -= ratio*a[k][j]
    return result


def mul(a,b):
    return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]


class EvolutionCertificates(unittest.TestCase):
    def test_iv_arbitrary_exact_blocks(self):
        steps = [[[2,3],[5,7]], [[-4,2],[8,9]], [[0,1],[1,0]]]
        a = defect_matrix(steps, 2)
        self.assertEqual(det(a), 1)
        self.assertEqual(det(mul(list(map(list,zip(*a))),a)), 1)

    def test_nilpotent_inverse_reconstructs_residuals(self):
        # Independent propagation formula applied to every residual basis vector.
        steps = [[[F(1, 2), 3], [2, -1]], [[0, 2], [F(3, 2), 1]]]
        l = defect_matrix(steps, 2)
        n = len(l)
        identity = [[int(i == j) for j in range(n)] for i in range(n)]
        inverse_columns = []
        for column in zip(*identity):
            state = list(column[:2])
            trajectory = list(state)
            for k, step in enumerate(steps):
                state = [column[2*(k+1)+i] + sum(a*b for a,b in zip(row, state))
                         for i,row in enumerate(step)]
                trajectory.extend(state)
            inverse_columns.append(trajectory)
        inverse = list(map(list, zip(*inverse_columns)))
        self.assertEqual(mul(l, inverse), identity)
        self.assertEqual(mul(inverse, l), identity)
        # Unit shift is nilpotent, independently confirming the finite inverse.
        shift = [[identity[i][j]-l[i][j] for j in range(n)] for i in range(n)]
        self.assertEqual(mul(shift, mul(shift, shift)), [[0]*n for _ in range(n)])

    def test_nonresonant_rotation_spectrum(self):
        # T=pi/2: M=r*J, so C^T C=(1+r^2)I, even at mu=0.
        for r in (F(1, 2), F(1), F(2)):
            monodromy = [[0, -r], [r, 0]]
            c = defect_matrix([monodromy], 2, periodic=True)
            self.assertEqual(mul(list(map(list, zip(*c))), c),
                             [[1+r*r, 0], [0, 1+r*r]])
            self.assertEqual(det(c), 1+r*r)

    def test_counterexample(self):
        a = defect_matrix([[[1]]], 1)
        self.assertEqual(a, [[1,0],[-1,1]])
        self.assertEqual(mul(list(map(list,zip(*a))),a), [[2,-1],[-1,1]])
        eigenvalues = [(3-math.sqrt(5))/2,(3+math.sqrt(5))/2]
        self.assertAlmostEqual(math.prod(eigenvalues), 1)
        self.assertAlmostEqual(positive_logdet(eigenvalues), 0)
        self.assertEqual(defect_matrix([],1), [[1]])

    def test_periodic_noncommuting_elimination(self):
        a,b,c = [[1,2],[0,1]], [[2,0],[3,1]], [[1,1],[1,0]]
        m = mul(c,mul(b,a))
        reduced = [[int(i==j)-m[i][j] for j in range(2)] for i in range(2)]
        self.assertEqual(det(defect_matrix([a,b,c],2,True)), det(reduced))

    def test_periodic_kernel_by_forward_substitution(self):
        a,b = [[1,1],[0,1]], [[1,-1],[0,1]]
        c = defect_matrix([a,b],2,True)
        for x in ([1,0],[0,1]):
            interior = [sum(v*w for v,w in zip(row,x)) for row in a]
            self.assertEqual([sum(v*w for v,w in zip(row,list(x)+interior)) for row in c], [0]*4)

    def test_rational_table(self):
        for r,delta,determinant in [(F(1,2),F(1,4),F(1,4)),(F(1),F(0),F(0)),(F(2),F(1),F(1))]:
            c = defect_matrix([[[r,0],[0,r]]],2,True)
            self.assertEqual(det(c), determinant)
            self.assertEqual(mul(list(map(list,zip(*c))),c), [[delta,0],[0,delta]])
        self.assertAlmostEqual(positive_logdet([.25,.25]), -4*math.log(2))

    def test_zero_and_scale_conventions(self):
        with self.assertRaises(ValueError): positive_logdet([0,0])
        self.assertEqual(positive_logdet([0,0],omit_zeros=True), 0)
        with self.assertRaises(ValueError): positive_logdet([-1])
        with self.assertRaises(ValueError): positive_logdet([float('nan')])
        self.assertLess(positive_logdet([1e-20,1]), 0)
        self.assertAlmostEqual(positive_logdet([2, 8]) - positive_logdet([1, 4]),
                               2*math.log(2))

    def test_signed_flow_and_square(self):
        # Finite-dimensional SF = endpoint negative inertia difference.
        start,end = [-1,1],[1,1]
        self.assertEqual(sum(x<0 for x in start)-sum(x<0 for x in end),1)
        self.assertEqual(sum(x*x<0 for x in start)-sum(x*x<0 for x in end),0)

    def test_variation_bound(self):
        # D(t)=exp(t), 0<=t<=1: true variation 2, gap 1.
        self.assertGreaterEqual(logdet_variation_bound(1,math.e-1,1),2)
        with self.assertRaises(ValueError): logdet_variation_bound(1,1,0)


if __name__ == '__main__': unittest.main()
