# Section 3: Corrected defect audit and Hopf normal-form illustration

Status (October 30, 2026): finite-dimensional algebra and a solvable ODE,
**not** a simulated or certified Navier–Stokes bifurcation.

## Initial-value defect: independent proof

The original equations define F(x)=Lx−b, b=(u₀,0,…,0), where
(Lx)₀=x₀ and (Lx)ᵢ=xᵢ−Sᵢxᵢ₋₁ for i=1,…,m−1.
For nonzero fixed u₀, F is affine. Its derivative L is the linear complex map.
For y=Lx, forward substitution gives x₀=y₀, xᵢ=yᵢ+Sᵢxᵢ₋₁. Therefore L is
invertible for every finite choice of Sᵢ. Block triangularity gives det L=1.
F=0 has a unique solution; it has no meaningful linear kernel.

With the same block metric in source and target, det(L*L)=1, all eigenvalues
are positive, Σlog λ=0 and λ_min≤1≤λ_max. Different source/target metrics add a
metric determinant ratio and must be declared; they do not create a kernel.
For fixed finite m and continuous finite Sᵢ on a compact parameter interval,
continuity gives a uniform positive gap. A small gap can indicate conditioning
or long-time amplification, not a finite kernel crossing. If ||Sᵢ||≤q<1,
write L=I−B with ||B||≤q to obtain σ_min(L)≥1−q. Contractivity is not necessary
for invertibility and was not established for the claimed NS discretization.

Minimal counterexample: m=1, H=R, L=[1], independently of any parameter.
A nontrivial scalar two-time example is L=[[1,0],[−1,1]], with
LᵀL=[[2,−1],[−1,1]], eigenvalues (3±√5)/2, determinant one, nullity zero,
and log determinant zero.

The old illustrative minima 2.4 and 1.1 are impossible for this operator. The
claimed real one-dimensional kernel, SF=+1, R=3520.3, Reynolds thresholds near
65/70 and receipt error/confidence/hash claims are withdrawn. Even the plausible
looking table rows were not measurements. If an eigenvalue of another operator
tends to zero while the others remain positive and bounded, Σlog λ tends to
**negative** infinity. At the singular point a primed sum excludes zeros and is
a different, possibly finite quantity. See Section 2 for regularization details.

## Omitted boundary conditions cannot be presumed

The initial row fixes x₀. No periodic condition was present. Dropping that row
leaves a rectangular recurrence with an H-dimensional kernel at every parameter.
A periodic alternative is (Cx)₀=x₀−S_m x_{m−1}, with the remaining recurrence
rows unchanged. Let M=S_m S_{m−1}…S₁. Eliminating interior variables by unit
triangular operations gives det C=det(I−M) and ker C≅ker(I−M).
For m=1, C=I−S₁ (both contributions occupy the same block).
Thus the periodic defect detects multiplier **one**, not every multiplier on
the unit circle. Autonomous periodic orbits already have a phase multiplier one;
a phase condition and often an unknown period are needed for continuation.

## Corrected exact example

Let J=[[0,−1],[1,0]], A(μ)=μI+J, T=2π. The exact equilibrium variational flow has
M(μ)=exp(TA)=exp(2πμ)I. Put r=exp(2πμ). Then C=(1−r)I and
Δ=CᵀC=(1−r)²I. Rational choices give exact matrices:

| r | μ | Δ spectrum | real nullity C | det C | R away from singularity |
|---|---|---|---:|---:|---|
| 1/2 | −log(2)/(2π) | 1/4, 1/4 | 0 | 1/4 | −4log(2) |
| 1 | 0 | 0, 0 | 2 | 0 | undefined; limiting value −∞ |
| 2 | log(2)/(2π) | 1, 1 | 0 | 1 | 0 |

At r=1 det' Δ=1 by the empty-product convention. The resonant real plane has
dimension two (one complex dimension). det C has no sign change and the positive
Laplacian has no signed SF. Choosing a nonresonant T can miss multiplier one at
μ=0 even when a Hopf bifurcation exists.

Now specify the nonlinear ODE z'=(μ+i)z−|z|²z. Polar coordinates give
ρ'=μρ−ρ³, θ'=1. For μ<0 the origin is asymptotically stable; for μ>0 it is
unstable with a stable cycle ρ=√μ of period 2π, radial derivative −2μ.
Eigenvalues μ±i cross with nonzero speed, and the cubic radial coefficient is −1.
These direct formulas prove this ODE's supercritical Hopf bifurcation. The
linear equation z'=(μ+i)z has the same equilibrium monodromy but no nonzero
periodic orbit for μ≠0. Monodromy singularity alone is insufficient.

For contrast B(t)=diag(t,1), t∈[−1,1], is self-adjoint with SF=+1 and real
nullity one at zero. Its square has eigenvalues t²,1 and loses the crossing sign.
This is a separate signed-crossing example, not a derived NS invariant.

## Reproduction and remaining physical proof obligations

Run `python -m unittest discover -s tests -v` and
`python examples/bifurcation_audit.py`. Rational algebra tests verify boundary
elimination and the replacement table; floating-point output independently
illustrates spectra/log sums. This is not interval-certified PDE evidence.

A Navier–Stokes certificate still needs domain and boundary conditions,
divergence-free basis and zero-mode handling, forcing and equilibrium branch,
convergence and validated spectral error bounds, an isolated conjugate pair with
transverse imaginary-axis crossing, exclusion of other critical spectrum, and
nonlinear nondegeneracy. A stochastic model requires a separately stated random
bifurcation notion. A fixed seed is not a deterministic Hopf theorem. None of
these physical proofs or data is supplied here. The former dimension count,
contractivity assertion and time-step specification were not validated evidence.
