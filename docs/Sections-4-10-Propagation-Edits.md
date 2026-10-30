# From operator certificates to geometric verification

A verification pipeline begins with a concrete question about the dynamics.
Given an initial state, how accurately has a trajectory been reconstructed?
Under periodic closure, which perturbations return to themselves? Across a
parameter path, how much can a spectral quantity change? The defect operator,
its boundary condition, and its metric turn each question into a matrix problem
with a specific certificate.

## Boundary choice and the size of the computation

Let the state space have dimension $d$ and the discretization contain $m$ state
blocks. The initial-value operator $L$ acts on $md$ coordinates and is invertible
by forward substitution. Its inverse quantifies how residuals accumulate through
the trajectory. Periodic closure gives an operator $C$ whose homogeneous solutions
are determined by their initial block. The return map $M$ gives the reduction

$$
\det C=\det(I-M),\qquad \ker C\cong\ker(I-M).
$$

Thus a periodic kernel computation on $md$ coordinates can be performed on the
$d\times d$ matrix $I-M$, followed by forward propagation of the resulting
initial vectors. The reduction preserves the determinant and kernel dimension.
Singular values retain the metrics and full propagation geometry; their bounds
are evaluated for the operator whose sensitivity is being certified.

For initial-value reconstruction, $\det L=1$ normalizes the product of singular
values. Sensitivity depends on their distribution: for example,
$\mathrm{diag}(K,K^{-1})$ has determinant one but amplifies some residuals by a
factor $K$ when $K>1$. This is why a receipt records both determinant data and a
lower singular-value bound.

## Residuals become state-error certificates

Suppose $Lx=b$ is the exact discrete equation, $\widehat x$ is a computed
trajectory, and $r=L\widehat x-b$. Since
$\widehat x-x=L^{-1}r$, a certified bound
$\sigma_{\min}(L)\ge s>0$ gives

$$
\|\widehat x-x\|_2\le\frac{\|r\|_2}{s}.
$$

A receipt containing $\|r\|_2\le\eta$ and $s$ therefore certifies the explicit
state-error ceiling $\eta/s$. The norm and gap use the declared metrics. When
floating-point arithmetic evaluates the residual, its rounding-error enclosure
is included in $\eta$. The estimate then measures error in the specified
discrete equation; transfer to a continuous model also uses its discretization
error bound.

The inputs and outputs are concrete:

| Certificate | Verified inputs | Output |
|---|---|---|
| Trajectory reconstruction | Operator and boundary data, residual ceiling $\eta$, singular gap $s$ | State error at most $\eta/s$ |
| Periodic kernel | Step maps and ordered monodromy product | Kernel of $I-M$, lifted through the recurrence |
| Approximate operator | Reference singular gap and operator-norm error | Invertibility margin and log-determinant error |
| Parameter variation | Uniform gap and integrated derivative bound | Change in log determinant over the path |

## Certification from an approximated operator

Let $\widehat D$ be an invertible $n\times n$ reference matrix with a verified
singular bound $\sigma_{\min}(\widehat D)\ge s>0$. Suppose the exact operator
$D$ satisfies the verified norm-error bound
$\|D-\widehat D\|_2\le\varepsilon<s$. The perturbation estimates in Section 2
give

$$
\sigma_{\min}(D)\ge s-\varepsilon>0,
$$

and, with $R(D)=\log\det(D^*D)$,

$$
|R(D)-R(\widehat D)|
\le -2n\log\left(1-\frac{\varepsilon}{s}\right).
$$

The positive margin $s-\varepsilon$ certifies invertibility throughout the
straight-line interpolation between the reference and exact operators. For real
matrices this preserves the determinant sign and its orientation class. It also
supplies a gap for residual estimates. The log-determinant estimate expresses
the same stability quantitatively: a small relative error
$\varepsilon/s$ gives a bound of order $2n\varepsilon/s$.

For a differentiable path of invertible $n\times n$ matrices with fixed metrics,
a uniform gap $\sigma_{\min}(D(t))\ge s>0$ instead gives

$$
|R(b)-R(a)|\le\frac{2n}{s}
\int_a^b\|D'(t)\|_2\,dt.
$$

Analytic estimates, interval enclosures, and validated residual methods can
supply these bounds. Lanczos and matrix-free estimators compute approximations
to the spectral quantities; the receipt carries the supporting enclosures that
turn those approximations into certificate inputs. Bounds evaluated in floating-point
arithmetic also include outward rounding in their final enclosures.

## Determinant-line transport

A smooth operator family supplies a family of determinant lines. A chosen
Hermitian metric and unitary connection determine geometric transport. For
$\nabla=d+A$, with curvature $F=dA$, an oriented spanning surface equipped with
a line trivialization satisfies

$$
\mathrm{Hol}(\partial\Sigma)
=\exp\left(-\int_\Sigma F\right).
$$

A geometric receipt accordingly specifies the parameter surface, connection,
orientation, and boundary path; open-path transport also specifies endpoint
frames. Signed spectral flow uses a self-adjoint family with declared crossing
and endpoint conventions. These data identify the geometric quantity being
computed alongside the magnitude and conditioning certificates above.

## Analytic and semantic extensions

The same certificate structure extends to broader settings once the following
inputs are established.

| Extension | Mathematical input required |
|---|---|
| Navier–Stokes Hopf bifurcation | Domain and boundary conditions, divergence-free state space, equilibrium branch, validated spectral convergence, transverse conjugate pair, nonlinear nondegeneracy |
| Infinite-dimensional determinants and anomalies | Operator domains, spectral hypotheses, zeta continuation regular at zero, normalization scale, geometric determinant-line construction |
| Stochastic geometry and entropy | Smooth parameter model, specified connection, proved relation to the chosen entropy functional |
| Witness cosheaves and semantic equivalence | Gluing laws, declared equivalence relation, soundness and completeness of the witnesses |
| Zero-knowledge verification | Explicit protocol with correctness and security proofs |
| Bifurcation forecasting | Quantitative relation between observed operator data and a future critical parameter |

The [example receipt](../receipts/schema.sample.json) records exact periodic
operator data and its norm-neighborhood certificate. Operator, spectral, and
transport receipts each have an explicit path from declared inputs to mathematical
conclusions.
