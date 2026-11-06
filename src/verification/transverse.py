"""Exact invariant restrictions; no orbit, section, or variational-ODE claims."""
import json
from .matrix import (matrix, MatrixEnclosure as Box, inverse_and_det, norm_upper,
                     periodic_return_enclosure, certify)
from .dynacomp_bridge import encode, sampled_step, replay_step


def identity(n):
    return matrix([[int(i == j) for j in range(n)] for i in range(n)])


def mul(a, b):
    a, b = matrix(a), matrix(b)
    if len(a[0]) != len(b):
        raise ValueError('product shape mismatch')
    return matrix([[sum(x*y for x, y in zip(row, col))
                    for col in zip(*b)] for row in a])


def frame(q, left):
    q, left = matrix(q), matrix(left)
    n, k = len(q), len(q[0])
    if not 0 < k < n or len(left) != k or len(left[0]) != n:
        raise ValueError('proper nonempty transverse frame required')
    if mul(left, q) != identity(k):
        raise ValueError('LQ must equal identity exactly')
    return q, left


def intertwiner(a, q, left):
    """Return B only after proving AQ=QB, not merely defining B=LAQ."""
    a = matrix(a)
    q, left = frame(q, left)
    if len(a) != len(q) or len(a[0]) != len(q):
        raise ValueError('ambient operator shape mismatch')
    aq = mul(a, q)
    b = mul(left, aq)
    if aq != mul(q, b):
        raise ValueError('chosen subspace is not exactly invariant')
    return b


def neutral_split(operators, q, neutral, *, generator):
    """Optional complete fixed neutral complement, verified over rationals."""
    neutral = matrix(neutral)
    n, k = len(q), len(q[0])
    if len(neutral) != n or len(neutral[0]) != n-k:
        raise ValueError('neutral complement shape mismatch')
    basis = matrix([list(r)+list(s) for r, s in zip(q, neutral)])
    inverse_and_det(basis)  # full direct sum, not a list of dependent vectors
    target = matrix([[0]*len(neutral[0]) for _ in range(n)]) if generator else neutral
    if any(mul(a, neutral) != target for a in operators):
        raise ValueError('complement is not pointwise neutral')


def _finish(q, left, boxes, neutral):
    cert = certify(periodic_return_enclosure(boxes))
    # ||QDz|| >= ||Dz||/||L|| and ||Qz|| <= ||Q|| ||z||.
    from fractions import Fraction
    ambient_gap = Fraction(cert['singular_gap_lower'])/(norm_upper(q)*norm_upper(left))
    return dict(certificate=cert, ambient_subspace_gap_lower=str(ambient_gap),
                metric='Euclidean in reduced coordinates; ambient bound separately',
                reduction='complete-fixed-neutral-complement' if neutral is not None
                          else 'invariant-subspace-only')


def exact_return_certificate(steps, q, left, *, neutral=None):
    """Exact rational discrete maps; same invariant frame at every step."""
    q, left = frame(q, left)
    steps = [matrix(a) for a in steps]
    if not steps:
        raise ValueError('at least one step required')
    reduced = [intertwiner(a, q, left) for a in steps]
    if neutral is not None:
        neutral = matrix(neutral)
        neutral_split(steps, q, neutral, generator=False)
    return dict(schema='quillen.transverse-exact/1', steps=[encode(a) for a in steps],
                Q=encode(q), L=encode(left), neutral=encode(neutral) if neutral is not None else None,
                reduced=[encode(a) for a in reduced],
                **_finish(q, left, [Box.point(a) for a in reduced], neutral))


def interval_restriction(box, q, left):
    """Sufficient universal invariance test for a supplied interval box.

    Containment of an intended physical map remains an external hypothesis.
    Zero lying in the residual box is insufficient; every endpoint must be zero.
    """
    q, left = frame(q, left)
    n = len(q)
    if len(box.lower) != n or len(box.lower[0]) != n:
        raise ValueError('ambient box shape mismatch')
    projection = mul(q, left)
    residual = Box.point([[int(i == j)-projection[i][j] for j in range(n)]
                          for i in range(n)]) @ box @ Box.point(q)
    if any(x != 0 for rows in (residual.lower, residual.upper) for row in rows for x in row):
        raise ValueError('box does not prove universal invariance')
    return Box.point(left) @ box @ Box.point(q)


def sampled_return_certificate(steps, q, left, *, neutral=None):
    """Source-bound DynaComp receipts + exact generator invariance proofs.

    The full and reduced exponential enclosures are independently produced using
    the merged bridge, then intersected. No interval residual proves invariance.
    """
    q, left = frame(q, left)
    steps = list(steps)
    if not steps:
        raise ValueError('at least one step required')
    operators, reduced, evidence, boxes = [], [], [], []
    for receipt in steps:
        full = replay_step(receipt)
        inv, _ = inverse_and_det(receipt['E'])
        a = mul(inv, receipt['A'])
        b = intertwiner(a, q, left)
        small, proof = sampled_step(b, receipt['h'])
        projected = Box.point(left) @ full @ Box.point(q)
        boxes.append(Box([[max(x,y) for x,y in zip(r,s)] for r,s in zip(small.lower,projected.lower)],
                         [[min(x,y) for x,y in zip(r,s)] for r,s in zip(small.upper,projected.upper)]))
        operators.append(a)
        reduced.append(encode(b))
        evidence.append(proof)
    if neutral is not None:
        neutral = matrix(neutral)
        neutral_split(operators, q, neutral, generator=True)
    return dict(schema='quillen.transverse-sampled/1', steps=steps, Q=encode(q), L=encode(left),
                neutral=encode(neutral) if neutral is not None else None,
                reduced_generators=reduced, reduced_steps=evidence,
                restricted_boxes=[dict(lower=encode(b.lower), upper=encode(b.upper)) for b in boxes],
                **_finish(q, left, boxes, neutral))


def replay(receipt):
    """Recompute all claims; unknown schemas and additional fields fail closed."""
    try:
        schemas = {'quillen.transverse-exact/1': exact_return_certificate,
                   'quillen.transverse-sampled/1': sampled_return_certificate}
        computed = schemas[receipt['schema']](receipt['steps'], receipt['Q'], receipt['L'],
                                              neutral=receipt['neutral'])
        return json.dumps(computed, sort_keys=True) == json.dumps(receipt, sort_keys=True)
    except (KeyError, TypeError, ValueError, ZeroDivisionError, IndexError):
        return False
