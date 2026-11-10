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


## Reading order and large-cat interference update

The cat tutorial now follows Exact Q Sampling as “Example: Cat States”, with
an explicit rejection-method introduction. The page URL and anchor remain
stable. The new cat-q-interference PNG/PDF shows α=3 analytic positive Q
functions on a shared log scale (colour floor 1e-9) and a Re z=0 line cut.
All curves share the constant reference exp(-9)/π=3.9282560692794874e-5;
unscaled density samples are in cat-q-interference-cut.tsv. There are no signed
weights or Wigner substitutions. The oscillatory factor period is π/3.

Validation: closed-form line cuts agree within 1e-12 in scaled units; exact
nodes are checked. An independent 80-level Fock calculation agrees within
1.275e-19 absolute Q error, with coherent-state omitted tail 4.237e-47. The
figures themselves have no Fock truncation. The even-cat central Q density is
7.856512018904454e-5, explaining its invisibility on the original linear scale.
Clean documentation build executes all seven tutorial blocks; HTTP checks
confirm the updated reading order and live interference section.

The rebuilt tutorial, new log-density/line-cut figure and reordered sidebar
were visually checked in the rendered local page at 1280 px width.


## Presentation audit: full plot versus central detail

The former large-cat lead image was a log-colour diagnostic with its central
curve divided by Qmix(0). This overemphasised the appearance of the very weak
central fringes relative to an ordinary full-Q plot. The replacement lead
image has common linear colour on [-5.5,5.5]² and raw-density central cuts.
The separate log diagnostic is preserved as cat-q-interference-log.png/pdf.
The central cut is correctly perpendicular to real ±α; its contrast scale is
1/(2α). Coordinate positions and Jacobians are documented explicitly.

Independent full-plane validation at 41×41 points (80-level Fock vectors) gives
maximum absolute Q error 2.498e-16. The analytic even-cat lobe height at z=α is
0.15915494551582188 and its origin density is 7.856512018904454e-5, only
0.049363920131% of that height. The central cut check retains 1.275e-19 absolute
agreement. The logarithmic plot and cut rescaling were a presentation issue,
not a change in the Q formula or sampler. No sampler source changed.
