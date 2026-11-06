# Exact linear propagators to periodic-return certificates

This optional adapter is staged on Quillen PR #2 at
`e9c8e7f14b0749c86f22a55719ae3b410f5cbf83`, which depends on PR #1.
It uses DynaComp main at `cd788c7ac695728e3f5b22056a0f99a2c1e1463c`.
DynaComp's default branch is separately `claude/new-repository-setup-yqzq7r`;
its Windows PR #1 is already merged there. No branches were integrated here.
No producer changes are needed, so its Explorer export contract is unaffected.

## Compatibility contract

`sampled_step(A,h,mass=E,stored=P)` accepts nonempty real square exact rational
matrices, positive exact rational h, and optionally a rational stored propagator.
Integer, Fraction and rational-string inputs are supported; floats and booleans
are rejected. States are columns in a single common coordinate system. E defaults
to identity. DynaComp performs exact descriptor elimination; singular E fails.
One dummy zero input makes its augmented dimension n+1, deliberately imposing
n<=11 and ||h E^-1 A||infinity<=8. No descriptor residual or numerical solve is
accepted. No uncertain coefficients, forcing-map return operator, changing basis,
nonlinear integration or time-varying segment is implied.

If P is omitted, a rational order-32 Taylor polynomial supplies a proposal only.
DynaComp independently computes its adaptive rational polynomial T_k and tail tau,
then returns an upward-rounded rational delta >= ||T_k-P||infinity+tau.
The exact endpoints P_ij +/- delta define Quillen's MatrixEnclosure. The row-norm
bound implies every entry bound since |X_ij| <= ||X||infinity. Thus the true
exp(h E^-1 A) belongs to this box. Stored coefficient discrepancy is included,
not added twice. The dummy B block is zero and does not influence the state map.
DynaComp's existing formats remain unchanged; bridge receipts use separate schemas.

For rho=||h[[E^-1 A,0],[0,0]]||infinity and k+2>rho, successive tail terms have
ratio at most rho/(k+2). Submultiplicativity and the geometric series give
 tau = rho^(k+1)/(k+1)! / (1-rho/(k+2)).
The producer seeks tau<=2^-100. In the scalar rho=1 case the order-32 majorant is
34/(33*33!)<2^-100. Exact zero rho is supported.

## Ordered return and determinant guarantees

`periodic_certificate(step_receipts)` replays source-bound producer evidence,
then calls the existing Quillen interval product and certificate APIs unchanged.
For chronologically ordered steps it encloses D=I-S_m ... S_1. Exact interval
multiplication is valid even for noncommuting matrices; independence of entries
is not assumed, and lost correlations merely enlarge the box.

Quillen computes an exact rational inverse and determinant of the box midpoint R.
Its rational norm majorants give s<=sigma_min(R), epsilon>=||D-R||2 and a positive
gap s-epsilon when epsilon<s. Its certificate also fixes determinant orientation
and bounds |log|det D|-log|det R|| by n epsilon/(s-epsilon). Hence
|det R| exp(-L) <= |det D| <= |det R| exp(L), with the certified sign.
These are determinant guarantees, not a new rational determinant-interval format.
A singular midpoint or insufficient margin raises rather than issuing a claim.

Receipts contain rational endpoints, original matrices, mass, interval, proposal,
DynaComp evidence and implementation hashes. Replay recomputes all of them and
Quillen's certificate. Hashes detect source drift, not producer authenticity or
formal correctness. The trusted base includes both projects' Python arithmetic,
source implementations and the stated mathematical arguments.

## Phase fixing is a separate problem

No neutral direction is removed automatically. If a periodic monodromy has a unit
multiplier, the full return defect may be singular and must fail certification.
The proposed [transverse certificate extension](Transverse-Return-Certificates.md)
checks exact generator/discrete-map invariance for a supplied rational basis Q
and left inverse L (LQ=I), then encloses the restricted defect I-L S Q using the
existing machinery. Optional complete fixed neutral complements are checked
separately, and coordinate/ambient metric bounds are distinguished. Merely
deleting one row/column does not prove this restriction.
For a Poincare section with normal n and orbit velocity f, the projection
I-f n^T/(n^T f) additionally requires a certified nonzero transversality margin,
validated f and n, return-time derivatives, and an established orbit. None is
supplied by this bridge. No nonlinear orbit existence is certified.

## Replay and validation

With sibling checkouts at the revisions above:

```sh
cd Quillen-Runtime
PYTHONPATH=../DynaComp/src python -m unittest discover -s tests -v
PYTHONPATH=../DynaComp/src python examples/dynacomp_return.py --replay receipts/dynacomp_return.json
PYTHONPATH=../DynaComp/src python examples/dynacomp_return.py > /tmp/bridge.json
diff receipts/dynacomp_return.json /tmp/bridge.json
cd ../DynaComp
PYTHONPATH=src python -m unittest discover -s tests
```

The independent scalar alternating-series test brackets exp(-1). Exact nilpotent
exponentials test noncommuting order without a numerical exponential oracle.
Retained examples include a coupled nonnormal stable segment and a noncommuting
two-segment return. Rejections exercise floats, booleans, shape errors, nonpositive
h, norm/dimension limits, singular E, tampering, a neutral map and excessive box
uncertainty. Optional bridge tests skip if DynaComp is absent; integration replay
must run with the pinned producer present.

Automatically certified: exact descriptor standardization, supported constant
segment exponential error including stored coefficients, step containment,
ordered-product containment, return-defect singular margin and determinant bounds.
External: exact source coefficients model the intended system; segments really
are constant and their durations exact; common coordinates and Euclidean metric
are appropriate; the selected segments form the desired period. Parameter
uncertainty, nonlinear orbit existence, variational linearization accuracy,
phase/section data and within-segment time variation remain external.

## Proposed publication operations (approval required)

Publish one new Quillen branch based on PR #2, containing only this adapter,
examples, tests, documentation and receipt. Open a dependent PR targeting
`codex/validated-matrix-certificates`. Merge neither existing PR nor DynaComp's
branches as part of this operation. After upstream integration, rebase/retarget
and rerun the pinned cross-project checks. No commits or remote writes have been
performed in this preparation.
