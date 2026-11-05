from fractions import Fraction as F
import importlib.util
import unittest
from src.verification.dynacomp_bridge import sampled_step, periodic_certificate, replay_periodic
from src.verification.matrix import periodic_return_enclosure

@unittest.skipUnless(importlib.util.find_spec('dynacomp'), 'optional DynaComp checkout required')
class BridgeTests(unittest.TestCase):
    def test_scalar(self):
        box,r=sampled_step([[-1]],1)
        term=total=F(1)
        for k in range(1,101):
            term*=F(-1,k); total+=term
        self.assertLessEqual(box.lower[0][0],total+term*F(-1,101))
        self.assertGreaterEqual(box.upper[0][0],total)
        c=periodic_certificate([r])
        self.assertGreater(F(c['certificate']['singular_gap_lower']),F(3,5))
        self.assertTrue(replay_periodic(c))
        c['steps'][0]['producer']['A_error_inf']='0'
        self.assertFalse(replay_periodic(c))
    def test_order(self):
        x,_=sampled_step([[0,1],[0,0]],1)
        y,_=sampled_step([[0,0],[1,0]],1)
        d=periodic_return_enclosure([x,y])
        for i,row in enumerate([[0,-1],[-1,-1]]):
            for j,v in enumerate(row):
                self.assertLessEqual(d.lower[i][j],v)
                self.assertGreaterEqual(d.upper[i][j],v)
        self.assertNotEqual(d.lower,periodic_return_enclosure([y,x]).lower)
    def test_descriptor(self):
        a,r=sampled_step([[-2]],1,mass=[[2]])
        b,_=sampled_step([[-1]],1,stored=[['1/3']])
        self.assertLessEqual(b.lower[0][0],a.lower[0][0])
        self.assertGreater(b.upper[0][0], F(367879,1000000))
        self.assertEqual(r['descriptor']['status'], 'exact')
    def test_rejections(self):
        cases=[([[0.1]],1,{}),([[True]],1,{}),([[1]],0,{}),([[9]],1,{}),
               ([[1,2]],1,{}),([[1]],1,{'mass':[[0]]}),
               ([[1]],1,{'stored':[[1,2]]}),
               ([[int(i==j) for j in range(12)] for i in range(12)],1,{})]
        for a,h,kw in cases:
            with self.subTest(a=a),self.assertRaises((ValueError,TypeError)):
                sampled_step(a,h,**kw)
        for a,kw in [([[0]],{}),([[-1]],{'stored':[[10]]})]:
            _,r=sampled_step(a,1,**kw)
            with self.assertRaises(ValueError): periodic_certificate([r])
