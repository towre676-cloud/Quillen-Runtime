# Computed matrix certificates

The matrix interface derives the inputs of Proposition 2.3 from a real square
matrix enclosure in fixed Euclidean metrics. It uses only Python's standard
library and exact rational arithmetic. It supports neither complex matrices nor
implicit changes of metric. Singular or insufficiently resolved inputs raise
`ValueError`; failure to certify is not a claim of singularity.

## Complete workflow

```python
from src.verification.matrix import MatrixEnclosure, certify, replay
box = MatrixEnclosure([['49/100', 0], [0, '49/100']],
                      [['51/100', 0], [0, '51/100']])
receipt = certify(box, candidate=[1, 0], rhs=['1/2', 0])
assert replay(receipt)
```

The reference defaults to the rational midpoint; an explicit reference need not
lie in the box. Entries accept integers, Fraction objects, or rational/decimal
strings. Floats are rejected: use `Fraction.from_float` explicitly for the exact
binary value, and supply measurement or integration uncertainty separately.

## Bounds and proof

For a rational reference A, exact pivoted elimination computes its inverse B and
determinant. For any matrix H define

$$U(H)=\operatorname{sqrtUpper}\left(\min\{\sum_{ij}H_{ij}^2,
\|H\|_1\|H\|_\infty\}\right).$$

Both quantities inside the minimum bound the squared Euclidean operator norm.
The upward square root is computed by integer square root on a dyadic grid of
spacing 2^-64, never by a floating point square root. Thus

$$s=1/U(B)\le\sigma_{\min}(A).$$

For interval endpoints L and V, put
$P_{ij}=\max(|L_{ij}-A_{ij}|,|V_{ij}-A_{ij}|)$.
Every enclosed D satisfies $|D-A|\le P$ entrywise and
$\|D-A\|_2\le\|P\|_2\le U(P)=\varepsilon$.
When $\varepsilon<s$, Proposition 2.3 yields the certified gap $s-\varepsilon$,
preserves the exact real determinant sign, and bounds

$$|\log\det(D^T D)-\log\det(A^T A)|
\le 2n\varepsilon/(s-\varepsilon).$$

This is a rational bound on the difference; it does not present floating point
log evaluations as certified absolute logarithmic determinant intervals.
The receipt gives det(A) exactly, so its reference value is $2\log|\det A|$.
For optional exact vectors x and b, the interface computes upper bounds for
$\|Ax-b\|_2$ and $\|x\|_2$ using the same arithmetic, then invokes the existing
solution-error theorem. The right hand side is fixed and exact in this version.

## Periodic variational systems

`periodic_return_enclosure(steps)` encloses $I-S_m\cdots S_1$ by exact outward
interval products, including noncommuting and nonnormal maps. Each step box must
already enclose the true variational propagator: a validated ODE integrator can
supply these endpoints, including its integration-error bounds. This interface
does not validate an ODE integration, orbit existence, or phase condition.
Autonomous periodic orbits generally have a neutral direction; a full return
map then cannot have a positive gap. Supply a justified transverse/phase-fixed
operator rather than silently discarding that mode. The certificate concerns
the reduced return defect, not a singular gap for the entire time-slice defect.

## Replay and limitations

```
python examples/validated_return_map.py
python examples/validated_return_map.py --replay receipts/validated_return_map.json
python -m unittest discover -s tests -v
```

Replay recomputes the inverse, determinant, norm bounds and theorem outputs from
the recorded rational data and requires exact agreement with every field. It is
an arithmetic check, not a signature or proof that the physical model lies in
the input box. Preserve the provenance of validated integration enclosures
separately. The implementation uses dense cubic elimination and conservative
norm bounds; interval wrapping and the fixed square-root grid may prevent
certification of a valid but poorly conditioned model. This release introduces
no eigenvalue-based gap assumptions. The nonnormal regression has eigenvalues
one and one but a singular gap below 1/100.
