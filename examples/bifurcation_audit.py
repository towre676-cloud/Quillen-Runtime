"""Standard-library floating-point illustration; exact checks live in tests."""
import json
import math
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from src.invariants.zeta import positive_logdet

rows=[]
for r in [.5,1.,2.]:
    lam=(1-r)**2
    rows.append({'r':r,'mu':math.log(r)/(2*math.pi),
                 'laplacian_eigenvalues':[lam,lam],
                 'real_nullity':2 if r==1 else 0,
                 'logdet':positive_logdet([lam,lam]) if lam else None,
                 'primed_logdet':positive_logdet([lam,lam],omit_zeros=True)})
print(json.dumps({'evidence':'floating_point_illustration_not_NS_simulation',
                  'counterexample_spectrum':[(3-math.sqrt(5))/2,(3+math.sqrt(5))/2],
                  'periodic_rows':rows},indent=2,allow_nan=False))
