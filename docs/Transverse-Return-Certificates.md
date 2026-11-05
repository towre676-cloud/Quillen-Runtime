# Certified transverse return maps for exact linear systems

This proposed extension is based on Quillen main
`4acfb7801266958cee50756ca411f6e0cd7dc043` and uses the existing DynaComp producer
at `cd788c7ac695728e3f5b22056a0f99a2c1e1463c`. It changes neither producer nor
the merged matrix-enclosure implementation. Its outputs are computer-assisted
linear-algebra certificates relative to those implementations and the arguments
below, not formal proof-assistant theorems or nonlinear orbit certificates.

## 1. Exact invariant-subspace theorem

Let Q be a rational n-by-k matrix, 0<k<n, and L a rational k-by-n matrix.
Require LQ=I exactly. This implies that Q has full column rank, identifies
coordinates on W=im(Q), and makes P=QL a projection onto W. It does not imply
orthogonality or dynamical invariance.

For an exact rational ambient operator A define B=LAQ. The additional identity

    AQ=QB

is equivalent to (I-QL)AQ=0 and establishes A-invariance of W. Defining B alone
does not establish that identity. By induction A^j Q=QB^j, including j=0.
Absolute convergence of the exponential series then gives

    exp(hA) Q = Q exp(hB).

For constant descriptor equations E x'=A x with invertible rational E, the
operator in this test is C=E^-1 A. Testing A rather than C is unsound. The
implementation uses the existing exact rational inverse, checks CQ=QB, and
replays the existing descriptor evidence; singular E fails.

For chronological segments j=1,...,m with a common invariant Q, set
S_j=exp(h_j C_j), R_j=exp(h_j B_j). No commutativity is required. Induction gives

    S_m ... S_1 Q = Q R_m ... R_1,
    D_W = I - R_m ... R_1 = I - L S_m ... S_1 Q.

This is the matrix of the restriction of I-S to W. The exact discrete-map API
uses the same proof with S_j Q=Q R_j directly, all over the rationals. Segment
invariance is sufficient, not necessary: a final product might preserve W even
when its factors do not. Such sampled systems are deliberately rejected here.
Changing frames, irrational invariant spaces and invariance known only for the
complete transcendental product are outside this API.

## 2. What “removing the neutral directions” proves

An invariant restriction by itself need not include every nonneutral direction.
Accordingly, an ordinary certificate says `invariant-subspace-only`.

An optional rational n-by-(n-k) matrix N supplies a stronger witness. The
implementation checks that [Q N] is invertible and that C_j N=0 for every
constant segment (or S_j N=N for exact discrete maps). Then S_j N=N and

    [Q N]^-1 S [Q N] = diag(R_m ... R_1, I).

If D_W is certified invertible, the full defect has kernel exactly im(N), rank
k, and precisely n-k unit multipliers counted algebraically. There is no hidden
unit eigenvalue or generalized unit block in the restricted part. This is a
complete fixed neutral complement, proved from the supplied exact operators.
The determinant reported is det(D_W); the full defect still has determinant
zero. A failure to find N is not proof that no invariant complement exists.

In particular, this option does not identify N with the tangent to a nonlinear
periodic orbit. All accepted constant segments annihilate the same N; many
legitimate periodic variational systems will not have this stronger property.

## 3. Restricted interval propagators and nonsingular defects

`sampled_return_certificate` replays every merged DynaComp step receipt. For each
reduced rational B_j it calls the same merged `sampled_step` implementation,
producing a second rigorous enclosure of exp(h_j B_j). It also projects the
ambient box with point matrices L and Q. By the invariance theorem, both boxes
contain the same reduced map, so their entrywise intersection does too. An empty
intersection fails. This approach preserves correlations proved symbolically
instead of trying to recover invariance from a noisy exponential box.

Existing `periodic_return_enclosure` multiplies the restricted boxes in the
correct order and encloses D_W. Existing `certify` supplies a rational reference
inverse and determinant, a positive lower singular margin, determinant sign and
log-determinant perturbation bound. Singular or insufficiently accurate defect
boxes fail. No eigenvalue-only test replaces a singular-value bound.

For reduced-coordinate gap g>0 the certificate proves

    ||D_W z||_2 >= g ||z||_2.

When Q is not orthonormal this is a coordinate Euclidean statement. Since
LQ=I, for x=Qz the ambient restriction obeys

    ||(I-S)x||_2 >= g/(||Q||_2 ||L||_2) ||x||_2.

The implementation outputs this conservative ambient bound using the existing
rational norm upper bounds. A change Q'=QT conjugates D_W, preserving determinant
and nonsingularity but generally changing its coordinate singular values.
No orthonormality is silently inferred from a rational basis.

The inherited determinant bound reports a rational upper bound for a logarithmic
error, not rational endpoints for the determinant. With reference determinant
d and output ell, the conservative enclosure
|d| exp(-ell) <= |det D_W| <= |d| exp(ell) and the certified sign are valid.
The inherited routine actually uses the factor-two bound for log det(D*D),
which is also a conservative bound for log|det D|; it is not sharpened here.

For externally supplied interval maps `interval_restriction` provides a separate
sufficient universal check: every endpoint of (I-QL)[S]Q must be exactly zero.
Then every matrix in the box preserves W, and L[S]Q is a valid restricted box.
Merely containing zero proves nothing. This test can reject valid correlated
families because interval arithmetic forgets correlations. Supplying a physical
map's enclosure remains an external mathematical hypothesis on this low-level
API; it is not promoted to source-bound DynaComp evidence.

## 4. The requested two-dimensional example

For C=diag(-1,0), h=1, Q=(1,0)^T, L=(1,0), N=(0,1)^T, every frame, invariance
and neutral-splitting identity is rational and exact. Thus

    S=diag(exp(-1),1),  D_W=1-exp(-1),  ker(I-S)=span(N).

The generated rational enclosure proves 1/2<D_W<2/3. Independently, the
alternating series gives 1/3<exp(-1)<1/2: the fourth partial sum is 3/8<1/2,
and the fifth is 11/30>1/3. The test uses the fortieth and forty-first partial
sums to cross-check the computed enclosure without a floating exponential.
The full defect is still singular and the original full-map certificate fails.

Two further retained examples exercise actual extensions. An oblique frame
Q=(2,2)^T, L=(1,-1/2) for C=[[-1,0],[-1,0]] has the same scalar reduced defect,
with coordinate gap approximately 0.6321205588 and conservative ambient gap
approximately 0.1998940722. A three-dimensional example uses two noncommuting
stable generators with transverse blocks

    B_1=[[-1,1],[0,-1]],  B_2=[[-1,0],[1,-1]]

and a third fixed neutral coordinate. For unit durations the reduced return is
exp(-2)[[1,1],[1,2]]. The computed rational defect certificate has gap greater
than 3/5 (approximately 0.6123097892 with the inherited norm estimate) and
positive determinant. Reversing the segments changes the certified matrix.
All decimal values here are illustrations of retained rational bounds.

## 5. Exact checks versus external assumptions

| Claim | Status in this implementation |
|---|---|
| Exact rational inputs, dimensions, E invertibility | Checked |
| LQ=I, CQ=QB or SQ=QR | Checked by exact equality |
| [Q N] invertible and fixed neutral complement | Checked when N is supplied |
| Exponential containment for supported constant segments | Replayed through the merged bridge |
| Restricted products, defect gap and determinant bounds | Computed through merged matrix machinery |
| General interval-family invariance | Sufficient universal zero-residual test only |
| Coefficients and durations model the intended physical system | External |
| Segments constitute the intended cycle in common coordinates | External |
| Invariance established by a separate analytic theorem | Not accepted as a boolean override; requires a future checked evidence interface |
| Irrational frames or uncertain coefficients preserving a correlated structure | Unsupported |
| Nonlinear periodic orbit and its tangent space | Not established |
| Validated variational equations along an orbit | Not established |
| Poincare return time and section derivative | Not established |

For a nonlinear autonomous orbit, a section normally requires a verified orbit,
normal n, nonzero n^T f, a return-time derivative and a projection such as
I-f n^T/(n^T f). A section tangent plane need not be invariant under fixed-time
monodromy. Compression by L and Q does not implement this projection or prove
the derivative of a first-return map. Neutral Jordan chains may prevent an
invariant complement altogether. The proposed APIs make none of these claims.

## 6. Reproduction and independent checking

With sibling Quillen and DynaComp checkouts at the pinned revisions:

```sh
PYTHONPATH=../DynaComp/src python -m unittest discover -s tests -v
PYTHONPATH=../DynaComp/src python examples/transverse_return.py > /tmp/transverse.json
diff receipts/transverse_return.json /tmp/transverse.json
PYTHONPATH=../DynaComp/src python examples/transverse_return.py --replay receipts/transverse_return.json
PYTHONPATH=../DynaComp/src python examples/check_transverse.py receipts/transverse_return.json
```

The separate checker never imports the new transverse module. It independently
implements the frame, intertwining and neutral-witness algebra, reconstructs the
restricted enclosures, and compares all claims. It intentionally reuses the
merged producer replay and matrix certificate machinery; it is independent of
the new reduction layer, not an independent proof of the entire numerical
stack. The trusted base still includes Python rational arithmetic, both merged
implementations, source-bound exponential evidence and the mathematical
arguments above. Hashes bind source versions, not formal correctness.

Rejections include row/column-deletion compressions with leakage, tiny exact
nonzero leakage, invalid left inverses, dependent neutral vectors, nonneutral
complements, zero-containing interval residuals, a singular reduced defect,
receipt tampering and unsupported fields. The tests also verify descriptor
counterexamples, oblique metrics, chronological ordering and independent
alternating-series bounds.
