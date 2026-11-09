# Cat-Q documentation and figure validation

Completed 10 November 2026 on `codex/exact-quantum-sampling`.

## Delivered

- `docs/src/cat-q.md`: visual tutorial covering ±α components, even/odd relative signs, a classical mixture, Q ordering, diagnostics, complex amplitudes, relative phases, loss and small odd cats.
- README, overview, state reference, support table and derivation links updated.
- `examples/cat_q.jl`: reproducible CairoMakie generator with independent analytic densities.
- PNG and PDF comparisons for α=1.2 and α=1.2 exp(iπ/6), plus the tab-separated validation record under `docs/src/assets`.

## Numerical checks

Julia 1.13.3 and CairoMakie 0.15.17; 300,000 samples per state, seeds 4242/4243. Six panels use a total of 1.8 million accepted samples. Each panel checks complex first and second moments, Q occupation, independently integrated density normalisation, probability in |Re z|<0.35 and |Im z|<0.35, and rejection acceptance. Six-standard-error tolerances apply to statistical checks. Central-region quadrature convergence is checked independently.

| α | State | Occupation error / SE | Central probability error / SE | Observed acceptance | Expected acceptance | Displayed mass |
|:--|:--|--:|--:|--:|--:|--:|
| 1.2 | Even cat | 1.74 | 1.89 | 0.528041 | 0.528067 | 0.999640 |
| 1.2 | Odd cat | 0.87 | 0.23 | 0.472611 | 0.471933 | 0.999667 |
| 1.2 | Classical mixture | 0.25 | 1.43 | 1.000000 | 1.000000 | 0.999657 |
| 1.2 exp(iπ/6) | Even cat | 0.00 | 0.25 | 0.527549 | 0.528067 | 0.999883 |
| 1.2 exp(iπ/6) | Odd cat | 1.27 | 2.15 | 0.472106 | 0.471933 | 0.999813 |
| 1.2 exp(iπ/6) | Classical mixture | 0.29 | 0.13 | 1.000000 | 1.000000 | 0.999837 |

The independent midpoint integrals over [-7,7]² differ from unit normalisation by less than 1e-12 for all six densities. Plotting uses 80×80 histogram bins over [-3.6,3.6]², analytic contours on a 401×401 grid, common axes and a shared density range [0,0.2]. Histograms divide by the complete sample count and retain the missing display-window mass; there is no spatial cutoff in sampling.

## Documentation and visual checks

- Clean Documenter HTML build passed with all 7 tutorial `@example` blocks executed; no warnings, broken cross references or errors.
- Package docstring doctests passed.
- Both standalone PNG comparisons inspected; main tutorial page inspected in headless Chrome at 1280×1100, including in-page figure/caption readability and equation rendering.
- PNG/PDF geometry is generated from the same Makie figures.
- `git diff --check` passed.
- No sampler implementation or package dependency changed; only documentation/example environments gained declarations for their imports. The existing full sampler test baseline remains recorded in `VALIDATION.md`.
- Local documentation build only; no deployment or publication.
