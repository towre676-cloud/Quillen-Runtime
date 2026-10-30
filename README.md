# Quillen-Runtime

A research scaffold for determinant-line ideas and computational receipts.
The October 2026 audit corrects the original, inconsistent Navier–Stokes example;
this repository does **not** currently certify a Navier–Stokes Hopf bifurcation
or equivalence of computational programs.

- [Finite-dimensional foundations and limits](docs/Section-2-Foundational-Definitions.md)
- [Defect-operator audit and corrected Hopf normal form](docs/Section-3-Navier-Stokes-Bifurcation.md)
- [Sections 4–10: surviving claims and proof obligations](docs/Sections-4-10-Propagation-Edits.md)
- [Illustrative receipt](receipts/schema.sample.json) (an example object, not a JSON Schema)

Run `python -m unittest discover -s tests -v` and
`python examples/bifurcation_audit.py` using Python 3.9+; no dependencies required.
Exact matrix checks and floating-point illustrations are explicitly distinguished.
No physical simulation, geometric holonomy implementation, or infinite-dimensional
zeta regularization is supplied.
