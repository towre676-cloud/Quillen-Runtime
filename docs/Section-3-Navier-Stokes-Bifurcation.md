# Section 3: Evolution defects, monodromy, and Hopf bifurcation

Boundary conditions determine which feature of a dynamical system a defect
operator measures. Initial-value defects encode reconstruction and amplification;
periodic defects encode return and resonance. This chapter gives exact formulas
for both, then connects periodic resonance to a nonlinear Hopf bifurcation whose
orbit, stability, and spectral data can be calculated in closed form.

## 3.1 Initial-value reconstruction and its determinant

Let $H$ be a $d$-dimensional real or complex inner-product space, and let
$S_1,\ldots,S_{m-1}:H\to H$ be linear evolution maps. On $H^{\oplus m}$ define

$$
(Lx)_0=x_0,\qquad (Lx)_i=x_i-S_i x_{i-1}\quad(1\le i<m).
$$

The prescribed initial state enters through the affine equation
$F(x)=Lx-b=0$, where $b=(u_0,0,\ldots,0)$. Its derivative is $L$.

**Proposition 3.1 (reconstruction and spectral normalization).** The operator
$L$ is invertible, and its inverse reconstructs the trajectory from its residuals:

$$
x_i=\sum_{j=0}^i\Phi_{i,j}y_j,\qquad
\Phi_{i,i}=I,\qquad \Phi_{i,j}=S_iS_{i-1}\cdots S_{j+1}\quad(j<i).
$$

With identical product metrics on source and target,

$$
\det L=1,\qquad \det(L^*L)=1,\qquad
\log\det(L^*L)=0,\qquad
0<\lambda_{\min}(L^*L)\le1\le\lambda_{\max}(L^*L).
$$

**Proof.** Write $L=I-B$, where $(Bx)_0=0$ and $(Bx)_i=S_ix_{i-1}$.
The block shift satisfies $B^m=0$, so

$$
L^{-1}=I+B+\cdots+B^{m-1}.
$$

Its blocks give the reconstruction formula. The matrix $L$ is block lower
triangular with identity diagonal, hence has determinant one. The remaining
identities follow from $\det(L^*L)=|\det L|^2$ and positivity: the product of
its $md$ eigenvalues is one. $\square$

The determinant normalizes total spectral volume, while the inverse measures
trajectory amplification. If $\|S_i\|_2\le q$ for $q\ge0$, then

$$
\|L^{-1}\|_2\le\Gamma_m(q):=\sum_{k=0}^{m-1}q^k,\qquad
\lambda_{\min}(L^*L)\ge\Gamma_m(q)^{-2}.
$$

This finite-horizon estimate includes the neutral case $q=1$, for which the gap
is at least $m^{-2}$. For $q<1$ it also gives the horizon-independent lower
bound $(1-q)^2$. On any compact parameter interval with continuous finite step
matrices, $L^{-1}$ remains bounded. Every member of the initial-value family is
joined to the identity through the invertible path $I-tB$, $0\le t\le1$.

For a computed trajectory $\widehat x$ with residual
$r=L\widehat x-b$, Proposition 3.1 therefore gives the a posteriori estimate

$$
\|\widehat x-x\|_2\le\Gamma_m(q)\|r\|_2.
$$

For example, $q\le1/2$ gives $\Gamma_m(q)<2$ for every finite horizon.
A residual norm at most $10^{-6}$ then bounds the trajectory error by
$2\cdot10^{-6}$, independently of the number of stored states. This turns the
inverse formula into a directly checkable accuracy criterion.

The smallest exact example is $m=d=1$, where $L=[1]$. The first example with
nontrivial propagation is

$$
L=\begin{pmatrix}1&0\\-1&1\end{pmatrix},\qquad
L^*L=\begin{pmatrix}2&-1\\-1&1\end{pmatrix},\qquad
\mathrm{spec}(L^*L)=\left\{\frac{3-\sqrt5}{2},\frac{3+\sqrt5}{2}\right\}.
$$

It has a positive spectral gap, zero log determinant, and reciprocal eigenvalues.
These identities provide an immediate consistency test for initial-value
receipts: a reported minimum above one or a nonzero log determinant contradicts
the specified operator and metric.

## 3.2 Periodic closure and monodromy reduction

For $m$ evolution maps $S_1,\ldots,S_m$, impose periodic closure on $m$ states:

$$
(Cx)_0=x_0-S_mx_{m-1},\qquad
(Cx)_i=x_i-S_ix_{i-1}\quad(1\le i<m).
$$

The monodromy map is the ordered product $M=S_mS_{m-1}\cdots S_1$.
For a single state, $C=I-S_1$.

**Proposition 3.2 (boundary reduction).** Periodic recurrence has the exact
finite-dimensional reduction

$$
\det C=\det(I-M),\qquad
\ker C\cong\ker(I-M).
$$

**Proof.** Use the residual coordinates $z_0=x_0$,
$z_i=x_i-S_ix_{i-1}$ for $i\ge1$. Their coordinate transformation is unit lower
triangular. Reconstructing $x_{m-1}$ in these coordinates gives

$$
(Cx)_0=(I-M)z_0-
\sum_{j=1}^{m-1}(S_mS_{m-1}\cdots S_{j+1})z_j,
\qquad (Cx)_i=z_i\quad(i\ge1).
$$

The transformed operator is block upper triangular with diagonal
$I-M,I,\ldots,I$. This proves the determinant identity. Its homogeneous
solutions are precisely the lifts

$$
v\longmapsto(v,S_1v,S_2S_1v,\ldots,S_{m-1}\cdots S_1v),
\qquad v\in\ker(I-M).
$$

The first component is the inverse of this kernel identification. $\square$

The boundary row controls the kernel: fixing the initial state gives a unique
trajectory; leaving it free gives $d$ homogeneous degrees of freedom; periodic
closure selects the fixed vectors of $M$. Consequently, a periodic defect
measures the multiplier-one eigenspace. For an autonomous periodic orbit the
tangent direction belongs to this eigenspace, motivating a phase condition when
solving for an orbit and its period.

The reduction moves determinant and nullspace calculations from an $md$-dimensional
trajectory space to the $d$-dimensional return map. Every reduced null vector
lifts explicitly to the full trajectory. This saves the construction and
factorization of the full block matrix when the target is a determinant or a
periodic mode. The reduction concerns determinants and kernels; conditioning
estimates also retain the norms of its coordinate transformations.

## 3.3 A complete resonance calculation

Identify $\mathbb C$ with $\mathbb R^2$, and consider the linear family

$$
A(\mu)=\mu I+J,\qquad
J=\begin{pmatrix}0&-1\\1&0\end{pmatrix},\qquad J^2=-I.
$$

Over a time interval $T>0$, its exact monodromy and one-state periodic defect are

$$
M_T(\mu)=e^{\mu T}(\cos T\,I+\sin T\,J),\qquad
C_T(\mu)=I-M_T(\mu).
$$

**Proposition 3.3 (resonance and multiplicity).** Put $r=e^{\mu T}$. Then

$$
C_T^\mathsf TC_T=\lambda(\mu,T)I,\qquad
\det C_T=\lambda(\mu,T),\qquad
\lambda(\mu,T)=1+r^2-2r\cos T.
$$

The kernel is two-dimensional over $\mathbb R$ exactly when
$\mu=0$ and $T\in2\pi\mathbb N$; otherwise $C_T$ is invertible.

**Proof.** Since $J^\mathsf T=-J$ and $J^2=-I$, the mixed terms cancel and

$$
\lambda(\mu,T)=(r-\cos T)^2+\sin^2T
               =|1-re^{iT}|^2.
$$

It vanishes exactly when $r=1$ and $e^{iT}=1$. At those parameters $C_T=0$.
The determinant follows by expanding the real $2\times2$ matrix. $\square$

For the natural rotation period $T=2\pi$, $C=(1-r)I$. Three rational choices
of $r$ give exact spectral data:

| $r$ | $\mu$ | Spectrum of $C^\mathsf TC$ | Real nullity | $\det C$ | Ordinary log determinant of $C^\mathsf TC$ |
|---|---|---|---:|---:|---|
| $1/2$ | $-\log 2/(2\pi)$ | $1/4,\;1/4$ | 0 | $1/4$ | $-4\log 2$ |
| $1$ | $0$ | $0,\;0$ | 2 | $0$ | Singular |
| $2$ | $\log 2/(2\pi)$ | $1,\;1$ | 0 | $1$ | $0$ |

Away from resonance,

$$
R(\mu)=4\log|1-e^{2\pi\mu}|
      =4\log|\mu|+4\log(2\pi)+O(\mu).
$$

Thus the two vanishing real singular directions give a logarithmic divergence
to $-\infty$. At the crossing the primed determinant omits both zero modes;
the empty product is one and its logarithm is zero. Section 2 gives this
convention and the corresponding metric and scale formulas.

## 3.4 The nonlinear Hopf branch and signed crossings

The nonlinear equation

$$
\dot z=(\mu+i)z-|z|^2z
$$

has the equilibrium linearization $A(\mu)$. Writing $z=\rho e^{i\theta}$ gives

$$
\dot\rho=\mu\rho-\rho^3,\qquad \dot\theta=1.
$$

The equilibrium eigenvalues $\mu\pm i$ cross the imaginary axis with unit
speed. For $\mu<0$ the origin is asymptotically stable. For $\mu>0$ it is
unstable and the exact periodic branch is

$$
z_\mu(t)=\sqrt\mu\,e^{i(t+\theta_0)}.
$$

Its period is $2\pi$. Linearizing the radial equation at $\rho=\sqrt\mu$ gives
$\dot\eta=-2\mu\eta$, so the cycle's Floquet multipliers are
$1$ in the phase direction and $e^{-4\pi\mu}$ in the radial direction.
The branch is orbitally asymptotically stable. These formulas establish the
supercritical Hopf bifurcation directly and identify both the resonant linear
plane at the equilibrium and the persistent phase direction along the cycle.

The linear equation $\dot z=(\mu+i)z$ shares the equilibrium monodromy, but its
nonzero amplitudes satisfy $\rho(t)=\rho(0)e^{\mu t}$. The cubic term supplies
the amplitude selection $\rho=\sqrt\mu$. This identifies the additional
nonlinear information used in the Hopf conclusion.

Signed spectral flow measures a different crossing. For the self-adjoint path

$$
B(t)=\mathrm{diag}(t,1),\qquad -1\le t\le1,
$$

the crossing form is positive on the one-dimensional kernel at zero, so
$\mathrm{SF}(B)=+1$. Its square has eigenvalues $t^2,1$ and spectral
flow zero. In the Hopf model, the Laplacian eigenvalues also meet zero from the
positive side: their singularity records resonance, while the equilibrium pair
$\mu\pm i$ records the direction of the stability change.

## 3.5 A quantitative return-map certificate

The rational case $r=1/2$ supplies the reference monodromy $\widehat M=I/2$,
two-dimensional reduced return defect $\widehat C=I/2$, and singular gap $s=1/2$. Suppose the monodromy
of interest satisfies the norm-error enclosure
$\|M-\widehat M\|_2\le1/100$. Since $C=I-M$, the same enclosure applies to
$C-\widehat C$. Proposition 2.3 yields

$$
\sigma_{\min}(C)\ge\frac{49}{100},\qquad
|R(C)-R(\widehat C)|\le-4\log(49/50)
\approx0.080811.
$$

The fully rational upper bound $4/49$ also encloses the log-determinant error.
A candidate reduced state with reference residual norm at most $10^{-6}$ and
state norm at most one consequently has error bounded by

$$
\|\widehat x-x\|_2
\le\frac{10^{-6}+1/100}{49/100}
=\frac{10001}{490000}\approx0.0204102.
$$

This calculation shows exactly how a return-map approximation becomes a
certificate: the monodromy error sets the invertibility margin and determinant
tolerance, while the residual and state norm set the solution tolerance. Every
real return defect in this norm ball also retains the positive determinant
orientation of the reference. In the
exact rational case the monodromy error is zero, reducing the state-error bound
to twice the residual norm. Near the resonant parameter, the formula identifies
the smaller operator error required to preserve a positive margin.

## 3.6 Reproduction

The standard-library implementation constructs both boundary operators. Exact
rational checks verify the determinant reduction, kernel lift, nilpotent inverse,
and spectral table. The numerical example evaluates the corresponding eigenvalues
and logarithms:

```bash
python -m unittest discover -s tests -v
python examples/hopf_spectra.py
```

The [operator receipt](../receipts/schema.sample.json) records the exact case
$r=1/2$. The [verification framework](Sections-4-10-Propagation-Edits.md) develops
these finite-dimensional formulas into derivative-and-gap certificates and
collects the analytic conditions for their transfer to Navier–Stokes.
