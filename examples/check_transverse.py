"""Independent transverse-layer checker; never imports verification.transverse.

Reuses the merged producer replay and matrix certificate implementation, whose
correctness remains in the trusted base. Exact frame/intertwining checks below
are separately implemented. Input is the bundle emitted by transverse_return.py.
"""
from fractions import Fraction as F
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.verification.matrix import (matrix, inverse_and_det, MatrixEnclosure as Box,
    periodic_return_enclosure, certify, norm_upper)
from src.verification.dynacomp_bridge import replay_step, sampled_step


def product(a, b):
    if len(a[0]) != len(b):
        raise ValueError('dimensions')
    return [[sum((a[i][t]*b[t][j] for t in range(len(b))), F(0))
             for j in range(len(b[0]))] for i in range(len(a))]


def equal(a, b):
    return matrix(a) == matrix(b)


def check(c):
    try:
        sampled = c['schema'] == 'quillen.transverse-sampled/1'
        if not sampled and c['schema'] != 'quillen.transverse-exact/1': return False
        keys = {'schema','steps','Q','L','neutral','certificate','ambient_subspace_gap_lower',
                'metric','reduction'} | ({'reduced_generators','reduced_steps','restricted_boxes'}
                                        if sampled else {'reduced'})
        if set(c) != keys: return False
        q,l=matrix(c['Q']),matrix(c['L'])
        n,k=len(q),len(q[0])
        if not 0<k<n or len(l)!=k or len(l[0])!=n: return False
        ident=[[F(i==j) for j in range(k)] for i in range(k)]
        if not equal(product(l,q),ident): return False
        if not c['steps']: return False
        bs,boxes,small_receipts=[],[],[]
        neutral=matrix(c['neutral']) if c['neutral'] is not None else None
        if neutral is not None:
            if len(neutral)!=n or len(neutral[0])!=n-k: return False
            inverse_and_det([list(q[i])+list(neutral[i]) for i in range(n)])
        for step in c['steps']:
            if sampled:
                full=replay_step(step)
                ei,_=inverse_and_det(step['E'])
                a=product(ei,matrix(step['A']))
            else:
                a=matrix(step)
            if len(a)!=n or len(a[0])!=n: return False
            aq=product(a,q); b=product(l,aq)
            if not equal(aq,product(q,b)): return False
            if neutral is not None:
                target=[[0]*(n-k) for _ in range(n)] if sampled else neutral
                if not equal(product(a,neutral),target): return False
            bs.append([[str(x) for x in row] for row in b])
            if sampled:
                small,sr=sampled_step(b,step['h']); small_receipts.append(sr)
                projected=Box.point(l) @ full @ Box.point(q)
                boxes.append(Box([[max(x,y) for x,y in zip(r,s)] for r,s in zip(small.lower,projected.lower)],
                                 [[min(x,y) for x,y in zip(r,s)] for r,s in zip(small.upper,projected.upper)]))
            else: boxes.append(Box.point(b))
        if bs != c['reduced_generators' if sampled else 'reduced']: return False
        if sampled:
            if small_receipts!=c['reduced_steps']: return False
            encoded=[dict(lower=[[str(x) for x in row] for row in b.lower],
                          upper=[[str(x) for x in row] for row in b.upper]) for b in boxes]
            if encoded!=c['restricted_boxes']: return False
        cert=certify(periodic_return_enclosure(boxes))
        if json.dumps(cert,sort_keys=True)!=json.dumps(c['certificate'],sort_keys=True): return False
        gap=F(cert['singular_gap_lower'])/(norm_upper(q)*norm_upper(l))
        if str(gap)!=c['ambient_subspace_gap_lower']: return False
        if c['metric']!='Euclidean in reduced coordinates; ambient bound separately': return False
        return c['reduction']==('complete-fixed-neutral-complement' if neutral is not None else 'invariant-subspace-only')
    except (KeyError,ValueError,TypeError,IndexError,ZeroDivisionError): return False


if __name__=='__main__':
    data=json.loads(Path(sys.argv[1]).read_text())
    if set(data)!={'diagonal','oblique','noncommuting'} or not all(check(c) for c in data.values()):
        raise SystemExit('independent transverse check rejected')
    print('independent transverse checks passed for all three examples')
