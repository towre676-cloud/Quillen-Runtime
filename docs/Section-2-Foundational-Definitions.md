# Section 2: Determinant lines, spectral variation, and transport

A finite evolution gives three complementary objects: a defect map, its positive
Laplacian, and a determinant line. The map records solvability, the Laplacian
measures conditioning, and a connection on the line defines geometric transport.
The results below specify these objects and the estimates connecting them.
Throughout, spaces are finite-dimensional over $\mathbb R$ or $\mathbb C$,
inner products are fixed unless stated otherwise, and $\|\cdot\|_2$ denotes the
induced operator norm.

## 2.1 The complex and its determinant line

For a linear map $D:V\to W$, consider the two-term complex
$$
0\longrightarrow V\xrightarrow{D}W\longrightarrow0,
\qquad H^0=\ker D,\qquad H^1=\operatorname{coker}D.
$$
Write $\det E=\Lambda^{\dim E}E$, with $\det\{0\}$ equal to the scalar field.
Our determinant-line convention is
$$
\operatorname{Det}H^\bullet(D)
 =\det(\ker D)\otimes\det(\operatorname{coker}D)^*.
$$
The exact sequence
$$
0\to\ker D\to V\to W\to\operatorname{coker}D\to0
$$
induces the determinant-functor identification
$$
\operatorname{Det}H^\bullet(D)\cong\det V\otimes(\det W)^*.
$$
This formulation also supplies a line through rank-changing families: the chain
spaces remain defined when the individual cohomology dimensions jump. For square
matrices, $\det D$ is naturally a section of the dual chain determinant line;
chosen bases make it an ordinary scalar. A metric and a connection are additional
structures. Complex Hermitian lines admit unitary phase transport; metric-compatible connections on real
one-dimensional bundles have orthogonal transport with values in
$\{\pm1\}$.

An affine defect $F(x)=Dx-b$ enters this construction through its linear part
$D$. The equation $F(x)=0$ describes solutions, while $\ker D$ describes
homogeneous variations. Section 3 applies this distinction to initial-value and
periodic boundary conditions.

## 2.2 Positive spectra and finite determinants

**Proposition 2.1 (paired spectra).** Set
$\Delta^0=D^*D$ and $\Delta^1=DD^*$. Their positive eigenvalues agree with
multiplicity, and
$$
\dim\ker\Delta^0=\dim V-\operatorname{rank}D,\qquad
\dim\ker\Delta^1=\dim W-\operatorname{rank}D.
$$
Consequently the index $\dim H^0-\dim H^1=\dim V-\dim W$ is independent of
$D$, including at rank changes.

*Proof.* In a singular-value decomposition $D=U\Sigma V^*$, the two
Laplacians are unitarily diagonalized by $V$ and $U$, respectively. Their
positive eigenvalues are the same squared singular values. The remaining entries
are zero, giving the stated dimensions. ∎

For a positive semidefinite operator $\Delta$, define
$$
\zeta_\Delta(s)=\sum_{\lambda_j>0}\lambda_j^{-s},\qquad
\det{}'\Delta=\prod_{\lambda_j>0}\lambda_j,\qquad
R(\Delta)=\log\det{}'\Delta=-\zeta_\Delta'(0).
$$
The sum is finite and entire in $s$. The empty product is one. At an invertible
operator the prime can be dropped. If a positive eigenvalue tends to zero while
the others stay positive and bounded, the ordinary logarithmic determinant tends
to $-\infty$; omitting the zero mode at the endpoint defines a different value.
For $c>0$,
$$
\log\det{}'(\Delta/c)=\log\det{}'\Delta
 -\operatorname{rank}(\Delta)\log c.
$$
Thus finite determinants require only a zero-mode convention and a declared
scale, without a cutoff or heat-kernel extrapolation.

For an invertible square map between $n$-dimensional spaces,
$$
\det(D^*D)=|\det D|^2
$$
in orthonormal bases. In arbitrary bases with Gram matrices $G_V,G_W$, the
formula becomes
$$
\det(D^*D)=|\det D|^2\frac{\det G_W}{\det G_V}.
$$
This separates the algebraic determinant from its metric normalization.

## 2.3 Quantitative variation of the logarithmic determinant

**Proposition 2.2 (variation and gap bounds).** Let $D(t)$, $a\le t\le b$,
be a $C^1$ path of invertible $n\times n$ matrices with fixed metrics. For
$R(t)=\log\det(D(t)^*D(t))$,
$$
R'(t)=2\operatorname{Re}\operatorname{tr}\bigl(D(t)^{-1}D'(t)\bigr).
$$
If $\sigma_{\min}(D(t))\ge s>0$ throughout the path, then
$$
|R(b)-R(a)|\le\frac{2n}{s}\int_a^b\|D'(t)\|_2\,dt.
$$
Alternatively, writing $\Delta(t)=D(t)^*D(t)$, a uniform bound
$\Delta(t)\ge gI$ gives
$$
|R(b)-R(a)|\le\frac{n}{g}\int_a^b\|\Delta'(t)\|_2\,dt.
$$

*Proof.* Jacobi's formula applied to $R=2\log|\det D|$ gives the derivative.
Use $|\operatorname{tr}A|\le n\|A\|_2$ and
$\|D^{-1}\|_2=1/\sigma_{\min}(D)$, then integrate. Applying the same argument
to $\log\det\Delta$ proves the second estimate. ∎

These bounds identify the data needed for a quantitative comparison: a uniform
gap and an integrated derivative bound. Gap estimates control sensitivity, while
the derivative integral controls the amount of variation along the path.

## 2.4 Projected connections and determinant transport

Let $P$ be a smooth constant-rank orthogonal projector onto a subbundle of a
trivial complex Hermitian bundle. A local orthonormal frame $Q$ satisfies
$Q^*Q=I$ and $P=QQ^*$. The projected connection and its induced determinant
connection have forms
$$
\mathcal A=Q^*dQ,\qquad A=\operatorname{tr}(Q^*dQ),\qquad
F=dA=\operatorname{tr}\bigl(P\,dP\wedge dP\bigr).
$$
The form $A$ is imaginary-valued. Under a unitary change of frame $Q\mapsto QU$,
it changes by $\operatorname{tr}(U^*dU)$, as a connection form should.
Differentiating $P^2=P$ gives $P(dP)P=0$, hence
$\operatorname{tr}(P\,dP)=0$. The frame expression above retains the phase data
that this identically vanishing trace cannot encode.

For a Hermitian line with $\nabla=d+A$, parallel transport is
$\exp(-\int A)$. If a closed curve bounds an oriented surface $\Sigma$ on
which a line trivialization is chosen, Stokes' theorem gives
$$
\operatorname{Hol}(\partial\Sigma)=\exp\left(-\int_\Sigma F\right),\qquad
\operatorname{dist}_{S^1}\bigl(\operatorname{Hol}(\partial\Sigma),1\bigr)
 \le\int_\Sigma|F|.
$$
Open-path scalar transport uses specified frames at both endpoints. This
curvature estimate controls phase around a boundary; Proposition 2.2 controls
logarithmic determinant magnitude along a path. For instance, $D(t)=[e^t]$
has $R(t)=2t$ even though every curvature two-form on its parameter interval
vanishes. A family indexed by random seeds becomes a curvature problem only after
a smooth parameter space and a connection have been supplied.

## 2.5 Signed crossings and spectral flow

For a continuous self-adjoint matrix path $B(t)$ with invertible endpoints,
choose the convention that a negative-to-positive crossing contributes positively.
Its spectral flow is
$$
\operatorname{SF}(B;[a,b])=n_-(B(a))-n_-(B(b)),
$$
where $n_-$ counts negative eigenvalues with multiplicity. For a $C^1$ path, at an isolated regular
crossing $t_0$, the contribution is the signature of the crossing form
$$
\Gamma_{t_0}(v,w)=\langle v,B'(t_0)w\rangle,
\qquad v,w\in\ker B(t_0).
$$
Each positive crossing-form direction removes one negative eigenvalue; each
negative direction adds one. Summing these changes proves the endpoint formula.
In particular every finite-dimensional closed loop has net spectral flow zero.

A positive Laplacian $D(t)^*D(t)$ has zero spectral flow whenever its endpoints
are invertible: an eigenvalue can touch zero without crossing from negative to
positive. For square $D$, the self-adjoint dilation
$$
\mathcal B(D)=\begin{pmatrix}0&D^*\\D&0\end{pmatrix}
$$
has paired eigenvalues $\pm\sigma_j(D)$ and also has zero net spectral flow
between invertible endpoints. Signed spectral flow, determinant-line holonomy and
complex determinant winding therefore encode different structures; relating them
requires a specified operator family and connection.

The [verification framework](Sections-4-10-Propagation-Edits.md) collects the analytic and geometric inputs for extending these results.
