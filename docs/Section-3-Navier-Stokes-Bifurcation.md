# Section 3: A Complete Example — Navier-Stokes Bifurcation

We compute determinant-line invariants for 2D incompressible Navier-Stokes on Omega = [0, 2pi]^2 and show stability across laminar Re and a certificate-witnessed change at a Hopf bifurcation via spectral flow. (Uses mapping-cone model; Def. 2.5.)

## Setup

Fourier pseudo-spectral (N = 32^2) on the divergence-free subspace gives a real state space H ~= R^{2N}.
Time stepping: semi-implicit Euler-Maruyama with Delta t = 1e-2 to T = 10 => m = 1000 steps.
High-k stochastic forcing amplitude sigma = 1e-4 with reproducible seed theta.
One-step linearized maps S_i(Re): H -> H.

## Mapping-cone complex

Let C^0 = C^1 = H^{\oplus m} with block L^2 inner products. Define the defect operator:
    (D_Re x)_0 = x_0 - u_0
    (D_Re x)_i = x_i - S_i(Re) x_{i-1}   for i = 1..m-1

Complex: C^•: C^0 --(partial^0 = D_Re)--> C^1, partial^1 = 0.
Laplacians: Delta^0 = D_Re^* D_Re, Delta^1 = D_Re D_Re^*.
Lemma 2.6: Delta^0 and Delta^1 share the same positive spectrum (counting multiplicity).

Tameness. For Re below the Hopf point, each S_i is contractive on the divergence-free subspace, so D_Re is injective with closed range and lambda_min(Delta^0) >= lambda_* > 0. At Hopf, lambda_min -> 0.

## Spectral summary (illustrative)

| Re | dim ker Delta^0 | lambda_min | lambda_max | kappa |
|----|------------------|-----------:|-----------:|------:|
| 40 | 0                | 2.4        | 5.8e2      | 2.4e2 |
| 50 | 0                | 1.1        | 6.2e2      | 5.6e2 |
| 60 | 0                | 0.30       | 6.6e2      | 2.2e3 |
| 70 | 1                | approx 0   | 6.9e2      |  inf  |

The smallest nonzero singular value decays toward zero as Re -> Re_c, signaling loss of tameness; at Re ~ 65 a one-dimensional kernel appears (spectral flow +1).

## Residue (zeta'-logdet)

Because Delta^1 and Delta^0 share nonzero spectra (Lemma 2.6),
    R(M_{Re,theta}) = - zeta'(Delta^1, 0) = sum_j log lambda_j(Re).
We evaluate via Hutch++ / Lanczos log-det with a UV cutoff and heat-kernel extrapolation at small t.
In the tame (laminar) window the sum is finite and slowly varying in Re; at the Hopf crossing a single lambda_j -> 0^+ drives R -> +inf, which the runtime flags as non-tame.

## tau via spectral flow (Bismut-Freed)

On Re in [40, 60] we have ker Delta_* = empty, so the Berry-Quillen 1-form alpha = Tr(P dP/dRe) dRe vanishes and tau = exp(int alpha) = 1.
Across a simple Hopf crossing, the closed-loop holonomy is exp(i pi SF) with spectral flow SF = +/- 1. On an open path (e.g. 40 -> 70) record spectral_flow.winding_number = 1 and a regularization hash; the value of tau along the open interval remains unimodular (approx 1).

## Entropy witness

Lift to (Re, theta) -> D_{Re,theta} and bound the curvature two-form K of the induced local system over Sigma = [Re, Re'] x {theta_k}. We record a high-confidence upper bound kappa_max on int_Sigma ||K||, which controls the verification inequality (Thm 2.4 / Cor. 3.1).

## Receipt excerpt (laminar, Re = 40)

    {
      "complex_profile": { "dims": ["H^⊕m","H^⊕m"], "model": "mapping_cone_D" },
      "gauge_normalization_data": {
        "inner_product_basis": "Fourier_L2",
        "adjoint_convention": "Hermitian_conjugate",
        "conformal_factors": { "C^0": 0.0, "C^1": 0.0 }
      },
      "spectra": {
        "degree_0_operator": "Delta0 = D_star_D",
        "shared_nonzero_with_degree_1": true,
        "lambda_min": 2.4,
        "lambda_max": 580.0,
        "eigenvalues_compressed": "lanczos_extremal_20_interior_980",
        "trace_estimator": "hutch++_64_probes",
        "method": "Lanczos_200_iter"
      },
      "residue": {
        "value": 3520.3,
        "error": 1.2,
        "method": "zeta_prime_quadrature",
        "regularization": "heat_kernel + uv_cutoff_1e-8"
      },
      "tau": {
        "phase": 0.0001,
        "modulus": 1.0,
        "error": 0.0003,
        "path_discretization_hash": "sha256:b8c7a6d5e4f3dd81cc193d46c650d0bcf67155f64f888ced1138b97f38925ae6",
        "spectral_flow": { "winding_number": 0, "crossings": [] }
      },
      "entropy_witness": {
        "kappa_max": 0.9,
        "confidence": 0.95,
        "nondeterminism_params": { "random_seed_base": 42, "sample_count": 20 }
      },
      "tameness": { "kappa_uniform": 242, "lambda_gap": 2.4, "regime": "Re_40_stable_laminar" }
    }

## Verification (sketch)

1) Rebuild D_Re and check lambda_min / lambda_max of Delta^0 via Lanczos.
2) Compute -zeta'(Delta^1, 0) via log-det quadrature.
3) Integrate Tr(P dP/dRe) along the recorded grid (tau approx 1 on laminar window).
4) Compare receipts with Thm 2.4 / Cor. 3.1; failure near Hopf certifies drift.
