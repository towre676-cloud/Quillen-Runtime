"""Optional DynaComp exact-linear sampling bridge. No numerical integration claims."""
from fractions import Fraction as F
import hashlib
import inspect
import json
from .matrix import matrix, rational, MatrixEnclosure, periodic_return_enclosure, certify


def encode(a):
    return [[str(x) for x in row] for row in a]


def source_digest():
    # Bind the actual producer implementation, not an unverified version label.
    import dynacomp.engineering as engineering
    import dynacomp.rational as algebra
    import dynacomp.machine as machine
    return {m.__name__: hashlib.sha256(inspect.getsource(m).encode()).hexdigest()
            for m in (engineering, algebra, machine)}


def sampled_step(a, h, *, mass=None, stored=None):
    """Enclose exp(h E^-1 A), for exact rational A,E,h and column states.

    One dummy zero input retains DynaComp's augmented sampling convention.
    This intentionally inherits its dimension/norm limits (n+1 <= 12).
    Floats, residual-only descriptor solves, and unsupported bounds fail closed.
    """
    from dynacomp.engineering import descriptor_transform, enclosure
    from dynacomp.rational import identity, mul
    a = matrix(a)
    n = len(a)
    if len(a[0]) != n or n+1 > 12:
        raise ValueError('square A with augmented dimension at most 12 required')
    h = rational(h)
    if h <= 0:
        raise ValueError('positive rational sampling interval required')
    e = matrix(mass) if mass is not None else identity(n)
    if len(e) != n or len(e[0]) != n:
        raise ValueError('mass shape mismatch')
    b = matrix([[0] for _ in range(n)])
    sa, sb, descriptor = descriptor_transform(a, b, e)
    if descriptor['status'] != 'exact':
        raise ValueError('exact descriptor elimination required')
    z = matrix([[h*x for x in row] for row in sa])
    if max(sum(abs(x) for x in row) for row in z) > 8:
        raise ValueError('augmented infinity norm exceeds 8')
    if stored is None:
        term = center = identity(n)
        for k in range(1, 33):
            term = matrix([[x/k for x in row] for row in mul(term, z)])
            center = matrix([[x+y for x,y in zip(r,s)] for r,s in zip(center,term)])
    else:
        center = matrix(stored)
        if len(center) != n or len(center[0]) != n:
            raise ValueError('stored propagator shape mismatch')
    evidence = enclosure(sa, sb, h, center, sb)
    if evidence['status'] != 'bounded':
        raise ValueError('DynaComp did not certify sampling')
    radius = F(evidence['A_error_inf'])
    box = MatrixEnclosure([[x-radius for x in row] for row in center],
                          [[x+radius for x in row] for row in center])
    receipt = dict(schema='quillen.dynacomp-step/1', A=encode(a), E=encode(e), h=str(h),
                   stored=encode(center), descriptor=descriptor, producer=evidence,
                   source_sha256=source_digest(), convention='column-state',
                   lower=encode(box.lower), upper=encode(box.upper))
    return box, receipt


def replay_step(receipt):
    box, computed = sampled_step(receipt['A'], receipt['h'], mass=receipt['E'],
                                 stored=receipt['stored'])
    if computed != receipt:
        raise ValueError('step receipt does not replay')
    return box


def periodic_certificate(steps):
    """Replay every producer receipt before independently certifying the defect."""
    steps = list(steps)
    result = certify(periodic_return_enclosure([replay_step(s) for s in steps]))
    return dict(schema='quillen.dynacomp-periodic/1', steps=steps, certificate=result)


def replay_periodic(receipt):
    try:
        return periodic_certificate(receipt['steps']) == receipt
    except (KeyError, TypeError, ValueError, ZeroDivisionError):
        return False
