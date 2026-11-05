# Quillen-Runtime

**Determinant lines, evolution operators, and quantitative spectral verification.**

Quillen-Runtime turns finite-dimensional evolution equations into reproducible
operator certificates. The boundary condition chooses the question: an
initial-value defect reconstructs a trajectory and measures sensitivity, while a
periodic defect tests which states return after one evolution cycle. Their
algebra is explicit:

$$
\det L=1,\qquad
\det C=\det(I-M),\qquad
\ker C\cong\ker(I-M).
$$

The periodic identity reduces a kernel problem on $m$ time slices of a
$d$-dimensional state space to the $d\times d$ return map $I-M$. For reconstruction,
a lower singular-value bound $s$ converts a computed residual into a state-error
bound:

$$
\|\widehat x-x\|_2\le
\frac{\|L\widehat x-b\|_2}{s},
\qquad Lx=b,\qquad \sigma_{\min}(L)\ge s>0.
$$

These quantities answer complementary questions. The determinant magnitude records the
product of singular values; the smallest singular value controls amplification
of errors. A determinant of one can therefore coexist with substantial trajectory
sensitivity. Verified matrix-error bounds extend the same estimates to an
approximated evolution operator and bound its log-determinant error.

A solvable Hopf normal form connects the operator formulas to nonlinear dynamics.
Its critical conjugate pair, emerging periodic orbit, and Floquet multipliers
are explicit. The accompanying treatment develops determinant-line conventions,
logarithmic spectral variation, geometric transport, and signed spectral flow.
Exact rational tests check the algebra, and a standard-library example reproduces
the spectral table.

| Chapter | Main results |
|---|---|
| [Determinant lines and spectral geometry](docs/Section-2-Foundational-Definitions.md) | Determinant conventions, paired spectra, residual and perturbation bounds, geometric transport, signed crossings |
| [Evolution defects and Hopf bifurcation](docs/Section-3-Navier-Stokes-Bifurcation.md) | Initial-value reconstruction, periodic dimension reduction, exact spectra, nonlinear orbit and stability |
| [From operator certificates to geometric verification](docs/Sections-4-10-Propagation-Edits.md) | Receipt inputs and outputs, numerical error certificates, and analytic extensions |

Run the verification suite and numerical example with Python 3.9 or later:

```bash
python -m unittest discover -s tests -v
python examples/hopf_spectra.py
```

Both commands use only the Python standard library. The
[example receipt](receipts/schema.sample.json) records exact matrix entries,
metrics, zero-mode conventions, and determinant data in a reproducible form.

## Validated matrix inputs

The [computed matrix certificate interface](docs/Validated-Matrix-Certificates.md)
derives spectral gaps, operator errors, residuals and state norms from exact
rational matrix enclosures. It produces replayable determinant and solution
certificates, including nonnormal periodic return maps from validated step boxes.

```bash
python examples/validated_return_map.py
python examples/validated_return_map.py --replay receipts/validated_return_map.json
```

## DynaComp sampled-step bridge

The [cross-project bridge](docs/DynaComp-Bridge.md) converts replayed bounded DynaComp descriptor/ZOH transformations into interval return maps. Exact noncommuting products feed the existing determinant and solution certificates. The pinned DynaComp producer is an optional dependency, tested explicitly by CI; unsupported/unbounded transformations are rejected.
