# From operator certificates to geometric verification

The identities in Sections 2 and 3 provide the algebraic core of a verification
pipeline. A dynamical model determines an evolution operator; boundary conditions
determine its defect complex; the complex supplies spectral quantities; a receipt
records the data needed to reproduce them. Geometric transport adds a connection
on the determinant line. Each stage has a distinct mathematical input and a
precise conclusion.

## Operator data and reproducible receipts

A receipt begins with the state space, field, ordered basis, source and target
metrics, evolution maps, and boundary condition. These data determine the adjoint,
the Laplacian, and the determinant convention. Initial-value and periodic defects
then have explicit algebraic certificates:

$$
\det L=1,\qquad
\det C=\det(I-M),\qquad
\ker C\cong\ker(I-M).
$$

Exact entries permit direct verification by rational elimination. For singular
operators the receipt records kernel dimension separately from the positive
spectrum and specifies whether its determinant omits zero modes. The
[example receipt](../receipts/schema.sample.json) uses this format for the
periodic normal form. Its matrix and spectral identities can be checked without
reconstructing an entire trajectory.

## Quantitative spectral verification

For a differentiable family of invertible $n\times n$ matrices with fixed metrics, write
$R(t)=\log\det(D(t)^*D(t))$. A uniform singular-value bound
$\sigma_{\min}(D(t))\ge s>0$ gives

$$
|R(b)-R(a)|\le \frac{2n}{s}
\int_a^b\|D'(t)\|_2\,dt.
$$

The certificate therefore has two numerical inputs: a lower bound valid along
the whole path and an upper bound on the derivative integral. Enclosures,
residual bounds, or analytic estimates can supply them. Lanczos quadrature and
matrix-free trace estimators offer computational routes to the spectral data;
their approximation errors belong in the receipt alongside the operator data.
A smallest Ritz value is an approximation to the smallest eigenvalue, whereas a
certified lower bound requires a supporting enclosure or estimate.

This structure separates exact identities, numerical approximations, and
validated error bounds. A verifier can check the relevant identity and then
propagate the declared bounds through the variation estimate.

## Determinant-line transport

A smooth operator family supplies a family of determinant lines. A chosen
Hermitian metric and unitary connection turn this algebraic object into a
geometric one. For $\nabla=d+A$, with curvature $F=dA$, transport around an
oriented spanning surface satisfies

$$
\operatorname{Hol}(\partial\Sigma)
 =\exp\left(-\int_\Sigma F\right).
$$

Accordingly, a geometric receipt specifies the parameter surface, connection,
orientation, and boundary path. Open-path transport additionally uses endpoint
frames. Spectral flow is recorded for a specified self-adjoint family, with its
crossing signs and endpoint conventions. These choices identify exactly which
geometric quantity a computation measures.

## Extension to analytic and semantic certificates

The finite-dimensional results establish operator and spectral certificates.
The following constructions extend that foundation to the broader programme.

| Extension | Mathematical input required |
|---|---|
| Navier–Stokes Hopf bifurcation | Domain and boundary conditions, divergence-free state space, equilibrium branch, validated spectral convergence, isolated transverse conjugate pair, and nonlinear nondegeneracy |
| Infinite-dimensional determinants and anomalies | Operator domains, spectral hypotheses, zeta continuation regular at zero, normalization scale, and geometric determinant-line construction |
| Stochastic geometry and entropy | Smooth parameter model, specified connection, and a proved relation between curvature and the chosen entropy functional |
| Witness cosheaves and program equivalence | Gluing laws, a declared semantic equivalence relation, and soundness and completeness results for the witnesses |
| Zero-knowledge verification | An explicit cryptographic protocol with correctness and security proofs |
| Bifurcation forecasting | A quantitative relation between observed operator data and the location of a future critical parameter |

These inputs determine the reach of each certificate: algebra identifies the
operator invariant, analysis transfers it to a dynamical system, and semantic or
cryptographic structure supplies the corresponding verification claim.
