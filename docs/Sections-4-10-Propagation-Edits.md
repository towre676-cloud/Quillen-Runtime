# Propagation Edits for Sections 4–10

Section 4 (Entropy as Curvature)
- Treat each run's cochain complex as the mapping-cone model (Definition 2.5).
- Families M_theta are families of D_{.,theta}. The curvature two-form K bounds variation through the parameter rectangle.
- The verification inequality (Theorem 2.4) uses an integral bound of ||K|| over Sigma.

Section 5 (Witness Cosheaf)
- At each x in the moduli space X, build C•(M) via the mapping-cone operator D.
- Receipts in W(U) are equivalence classes of (tau, R) modulo the verification inequality.

Section 6 (Anomaly Theory)
- Determinant gerbe over mapping-cone complexes.
- Delta^0 = D*D and Delta^1 = D D* share the same nonzero spectrum (Lemma 2.6).
- The chosen zeta regularization is part of the trivialization data.

Section 7 (Operator Theory)
- Delta^0 encodes the discrete evolution energy derived from D.
- tau is computed via the Berry-Quillen connection on the determinant line.
- R is computed via zeta on Delta^0 or Delta^1 (equivalent by Lemma 2.6).

Section 8 (Implementation)
- Under tameness (uniform spectral gap and condition number bounds), (tau, R) is Lipschitz in the parameter.
- Computation to tolerance epsilon is done by Lanczos-based spectral quadrature.
- Verification is equivalent to the inequalities of Theorem 2.4.

Section 9 (Zero-Knowledge)
- Matrix-free Delta^0 log-det; bound entropy using Delta^1 without revealing D.

Section 10 (Predictive Holonomy)
- Curvature norms grow near spectral-flow crossings (see Section 3).
- Use Corollary 3.1 to forecast the crossing before reaching it.

Conclusion add-on
- The invariance theorems, mapping-cone construction, and the verification inequality reduce equivalence checking to finite numerical bounds rather than whole-program re-execution.
