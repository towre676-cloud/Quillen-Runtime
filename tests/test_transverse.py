import copy
import importlib.util
import itertools
import json
import unittest
from fractions import Fraction as F
from src.verification.matrix import MatrixEnclosure as Box, certify, periodic_return_enclosure
from src.verification.dynacomp_bridge import sampled_step
from src.verification.transverse import (frame, intertwiner, exact_return_certificate,
    sampled_return_certificate, interval_restriction, replay)

Q = [[1], [0]]
L = [[1, 0]]
N = [[0], [1]]


class ExactTransverse(unittest.TestCase):
    def test_exact_neutral_split(self):
        c = exact_return_certificate([[['1/3', 0], [0, 1]]], Q, L, neutral=N)
        self.assertEqual(F(c['certificate']['singular_gap_lower']), F(2, 3))
        self.assertEqual(c['reduction'], 'complete-fixed-neutral-complement')
        self.assertTrue(replay(json.loads(json.dumps(c))))

    def test_oblique_frame_and_metric(self):
        # Q is not unit length; L is a valid oblique left inverse.
        c = exact_return_certificate([[['1/2', 0], ['-1/2', 1]]], [[2], [2]], [[1, '-1/2']], neutral=N)
        self.assertEqual(c['reduced'], [[['1/2']]])
        self.assertLess(F(c['ambient_subspace_gap_lower']), F(1, 2))

    def test_compression_without_invariance_rejected(self):
        # Deleting the second row/column suggests gap 1/2 but e1 is not invariant.
        with self.assertRaises(ValueError):
            exact_return_certificate([[['1/2', 0], [1, 1]]], Q, L)

    def test_noncommuting_order(self):
        a = [[1, 1, 0], [0, 1, 0], [0, 0, 1]]
        b = [[2, 0, 0], [2, 2, 0], [0, 0, 1]]
        q, l, n = [[1,0],[0,1],[0,0]], [[1,0,0],[0,1,0]], [[0],[0],[1]]
        x = exact_return_certificate([a,b],q,l,neutral=n)
        y = exact_return_certificate([b,a],q,l,neutral=n)
        self.assertEqual(x['certificate']['lower'], [['-1','-2'],['-2','-3']])
        self.assertNotEqual(x['certificate']['lower'], y['certificate']['lower'])

    def test_invalid_frames(self):
        for q,l in [([[0],[0]],L), (Q,[[2,0]]), ([[1,0],[0,1]],[[1,0],[0,1]]),
                    (Q,[[1]]), ([[1.0],[0]],L), ([[True],[0]],L)]:
            with self.subTest(q=q), self.assertRaises(ValueError): frame(q,l)

    def test_invalid_neutral_complements(self):
        for n in (Q, [[1],[1]], [[0,0],[1,1]]):
            with self.subTest(n=n), self.assertRaises(ValueError):
                exact_return_certificate([[['1/2',0],[0,1]]],Q,L,neutral=n)

    def test_invariant_only_does_not_claim_neutral(self):
        c=exact_return_certificate([[['1/2',1],[0,2]]],Q,L)
        self.assertEqual(c['reduction'],'invariant-subspace-only')

    def test_singular_restriction_and_empty_fail(self):
        for maps in ([], [[[1,0],[0,1]]]):
            with self.assertRaises(ValueError): exact_return_certificate(maps,Q,L)

    def test_interval_universal_invariance(self):
        b=Box([['1/3',-2],[0,1]],[['1/2',3],[0,2]])
        r=interval_restriction(b,Q,L)
        self.assertEqual(r.lower,((F(1,3),),))
        self.assertEqual(r.upper,((F(1,2),),))
        # Zero-containing residual does not prove exact invariance.
        with self.assertRaises(ValueError):
            interval_restriction(Box([['1/3',0],[-1,1]],[['1/2',0],[1,1]]),Q,L)

    def test_interval_all_vertices(self):
        for lo,hi in [(F(1,4),F(1,2)),(F(-1,2),F(1,3))]:
            r=interval_restriction(Box([[lo,-1],[0,1]],[[hi,1],[0,1]]),Q,L)
            for x,y in itertools.product((lo,hi),(-1,1)):
                self.assertLessEqual(r.lower[0][0],x)
                self.assertGreaterEqual(r.upper[0][0],x)

    def test_tampering_all_fields(self):
        c=exact_return_certificate([[['1/2',0],[0,1]]],Q,L,neutral=N)
        for key in c:
            d=copy.deepcopy(c); d[key]=None
            self.assertFalse(replay(d),key)
        c['external_invariance_assumed']=True
        self.assertFalse(replay(c))
        for obj in (None, [], {}, {'schema':'unknown'}): self.assertFalse(replay(obj))


@unittest.skipUnless(importlib.util.find_spec('dynacomp'),'pinned DynaComp required')
class SampledTransverse(unittest.TestCase):
    def test_requested_example_independent_alternating_series(self):
        full,step=sampled_step([[-1,0],[0,0]],1)
        with self.assertRaises(ValueError): certify(periodic_return_enclosure([full]))
        c=sampled_return_certificate([step],Q,L,neutral=N)
        lo=F(c['certificate']['lower'][0][0]); hi=F(c['certificate']['upper'][0][0])
        self.assertGreater(lo,F(1,2)); self.assertLess(hi,F(2,3))
        # Independent exact alternating-series enclosure, no exp oracle.
        term=total=F(1)
        for k in range(1,41): term*=F(-1,k); total+=term
        lower_exp=total+term*F(-1,41)
        self.assertLessEqual(lo,1-total)
        self.assertGreaterEqual(hi,1-lower_exp)
        self.assertTrue(replay(json.loads(json.dumps(c))))

    def test_descriptor_standardization_is_essential(self):
        # A itself does not preserve e1; E^-1 A does.
        _,s=sampled_step([[-2,0],[-1,0]],1,mass=[[2,0],[1,1]])
        c=sampled_return_certificate([s],Q,L,neutral=N)
        self.assertEqual(c['reduced_generators'],[[['-1']]])
        # Conversely A preserves e1 while E^-1 A does not.
        _,s=sampled_step([[-1,0],[0,0]],1,mass=[[1,0],[1,1]])
        with self.assertRaises(ValueError): sampled_return_certificate([s],Q,L)

    def test_noncommuting_generators(self):
        a=[[-1,1,0],[0,-1,0],[0,0,0]]
        b=[[-1,0,0],[1,-1,0],[0,0,0]]
        q,l,n=[[1,0],[0,1],[0,0]],[[1,0,0],[0,1,0]],[[0],[0],[1]]
        steps=[sampled_step(x,1)[1] for x in (a,b)]
        c=sampled_return_certificate(steps,q,l,neutral=n)
        d=sampled_return_certificate(steps[::-1],q,l,neutral=n)
        self.assertTrue(replay(c))
        self.assertGreater(F(c['certificate']['singular_gap_lower']),F(1,2))
        self.assertNotEqual(c['certificate']['lower'],d['certificate']['lower'])

    def test_oblique_generator(self):
        _,s=sampled_step([[-1,0],[-1,0]],1)
        c=sampled_return_certificate([s],[[2],[2]],[[1,'-1/2']],neutral=N)
        self.assertTrue(replay(c))
        self.assertEqual(c['reduced_generators'],[[['-1']]])

    def test_small_nonzero_leakage_fails(self):
        _,s=sampled_step([[-1,0],['1/100000000000000000000',0]],1)
        with self.assertRaises(ValueError): sampled_return_certificate([s],Q,L)

    def test_neutral_claim_and_singular_restriction_fail(self):
        _,s=sampled_step([[-1,0],[0,1]],1)
        with self.assertRaises(ValueError): sampled_return_certificate([s],Q,L,neutral=N)
        _,s=sampled_step([[0,0],[0,0]],1)
        with self.assertRaises(ValueError): sampled_return_certificate([s],Q,L,neutral=N)

    def test_tampering_all_fields(self):
        _,s=sampled_step([[-1,0],[0,0]],1)
        c=sampled_return_certificate([s],Q,L,neutral=N)
        for key in c:
            d=copy.deepcopy(c); d[key]=None
            self.assertFalse(replay(d),key)
        d=copy.deepcopy(c); d['steps'][0]['A'][1][0]='1/1000'
        self.assertFalse(replay(d))
        d=copy.deepcopy(c); d['certificate']['singular_gap_lower']='1'
        self.assertFalse(replay(d))

    def test_separate_checker_and_tampering(self):
        from examples.check_transverse import check
        from examples.transverse_return import examples
        for c in examples().values():
            self.assertTrue(check(c))
            for key in c:
                d=copy.deepcopy(c); d[key]=None
                self.assertFalse(check(d),key)
        c=exact_return_certificate([[['1/2',0],[0,1]]],Q,L,neutral=N)
        self.assertTrue(check(c))
