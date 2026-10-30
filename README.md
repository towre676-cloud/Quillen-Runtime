# Quillen-Runtime

**Determinant lines, evolution operators, and quantitative spectral verification.**

Quillen-Runtime develops a finite-dimensional framework for extracting spectral
certificates from discrete dynamics. Its central result identifies how boundary
conditions determine the invariant: an initial-value defect has determinant one,
while periodic closure reduces its determinant and kernel to the monodromy map.

$$
\det L=1,\qquad
\det C=\det(I-M),\qquad
\ker C\cong\ker(I-M).
$$

The treatment connects these identities to determinant lines, log determinants,
conditioning, and signed spectral flow. A fully solvable Hopf normal form supplies
the nonlinear dynamics: the critical conjugate pair, the emerging periodic orbit,
and its stability are all explicit. Exact rational matrix checks accompany the
algebra, and a standard-library example reproduces the spectral table.

| Chapter | Main results |
|---|---|
| [Determinant lines and spectral geometry](docs/Section-2-Foundational-Definitions.md) | Cohomological determinant conventions, positive-spectrum identity, log-determinant variation, connection and spectral-flow formulas |
| [Evolution defects and Hopf bifurcation](docs/Section-3-Navier-Stokes-Bifurcation.md) | Initial-value invertibility, periodic monodromy reduction, exact spectra, nonlinear normal form |
| [From operator certificates to geometric verification](docs/Sections-4-10-Propagation-Edits.md) | Receipt structure, certified numerical bounds, geometric transport, and analytic extensions |

Run the verification suite and numerical example with Python 3.9 or later:

```bash
python -m unittest discover -s tests -v
python examples/hopf_spectra.py
```

Both commands use only the Python standard library. The
[example receipt](receipts/schema.sample.json) records exact matrix entries,
metrics, zero-mode conventions, and determinant data in a reproducible form.
