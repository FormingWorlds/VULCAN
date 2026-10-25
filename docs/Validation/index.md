# Validation

## What this page checks

This tab lists the tests that compare VULCAN with an independent result: a published benchmark, an analytical limit or a second implementation. Validation means that a test runs the code on a case with a known right answer and checks the result within a stated tolerance. A wrong result in one of these parts of the code would change the chemistry that VULCAN computes for the atmosphere. Each case page says why the check exists, what it tests, how we know the right answer and why it matters. A box "Test details" at the end of each case gives its reference, its tolerance, its test id and the command that runs it.

## Cases by process

| Process | Case | Test id | Reference | Status |
|---|---|---|---|---|
| [Stellar irradiation](stellar_irradiation.md) | Top-of-atmosphere flux | `tests/test_integration_quick.py::test_top_of_atmosphere_flux_is_diluted_stellar_spectrum` | stellar input file diluted by (R_star / a)^2, interpolated with numpy | integration tier; plot pair |
| Chemical kinetics | none yet | | | no published or analytical benchmark of the chemistry is tested yet |

## Numerical checks

Two tests compare the linear solvers with `numpy.linalg.solve` on the same matrices, to `rtol=1e-12` (the Thomas test also with `atol=1e-12`, the band test with `atol=0`). They check the numerics of the solver, not the physics, and are not marked as validation cases:

- `tests/test_op.py::test_thomas_vec_matches_dense_solve`: the Thomas algorithm for a tridiagonal system with several right-hand sides.
- `tests/test_op.py::test_band_solve_matches_dense_solve`: the band LU factorisation of the Ros2 solver, with chemistry off.

## Plots written by the tests

A case with a plot pair shows the tested values against their reference in a light and a dark version. The test writes the pair with the same numbers it asserts, through the `validation_plot` fixture of proteus-mpl. This command writes all pairs to `docs/assets/validation/`:

```console
pytest -m reference_pinned tests --validation-plots
```

Without `--validation-plots` the tests assert and write nothing. The command makes every plot again; commit only the pair of the case you changed, because the bytes match only for the same matplotlib and Pillow versions. The integration cases need the FastChem executable (`make` in `fastchem_vulcan/`).
