# Final validation

10 November 2026. Starting revision b180b50f0eddfed02d2f19728c6cbb5d15865b36.
Both runtimes are existing juliaup-managed channels; no launcher, depot or channel was replaced.

## Julia 1.10.12

Full `Pkg.test("PhaseSpaceTools")`: 1572 assertions passed across 32 top-level testsets.

- Scientific regressions: 104/104
- Aqua: 11/11
- Arg tests: 69/69
- Fock W warning: 2/2
- Noises: 6/6
- Squeezed Q: 3/3
- Crescent Q: 3/3
- Coherent W: 3/3
- Squeezed2 Q: 4/4
- Coherent +P: 4/4
- Fock +W small n: 3/3
- Fock +W asymptotic: 3/3
- Thermal P,Q,W: 3/3
- Fock W: 4/4
- Fock +W: 4/4
- Fock +W (n≫1): 4/4
- Fock +P: 6/6
- Squeezed +P: 2/2
- Squeezed2 +P: 4/4
- Squeezed W: 2/2
- Bogoliubov W: 3/3
- Bogoliubov Q: 3/3
- Bogoliubov +P: 3/3
- RNG, layout and validation: 102/102
- Exact Q families: 104/104
- Finite superpositions and cats: 189/189
- Independent odd-cat density and envelope: 90/90
- Passive mixing and collective Fock states: 349/349
- Loss and phase diffusion: 99/99
- Gaussian covariances and quadrature conventions: 169/169
- Photon-added centred Gaussian: 77/77
- Exact doubled positive-W: 140/140

## Julia 1.13.3

Full `Pkg.test("PhaseSpaceTools")`: 1572 assertions passed across 32 top-level testsets.

- Scientific regressions: 104/104
- Aqua: 11/11
- Arg tests: 69/69
- Fock W warning: 2/2
- Noises: 6/6
- Squeezed Q: 3/3
- Crescent Q: 3/3
- Coherent W: 3/3
- Squeezed2 Q: 4/4
- Coherent +P: 4/4
- Fock +W small n: 3/3
- Fock +W asymptotic: 3/3
- Thermal P,Q,W: 3/3
- Fock W: 4/4
- Fock +W: 4/4
- Fock +W (n≫1): 4/4
- Fock +P: 6/6
- Squeezed +P: 2/2
- Squeezed2 +P: 4/4
- Squeezed W: 2/2
- Bogoliubov W: 3/3
- Bogoliubov Q: 3/3
- Bogoliubov +P: 3/3
- RNG, layout and validation: 102/102
- Exact Q families: 104/104
- Finite superpositions and cats: 189/189
- Independent odd-cat density and envelope: 90/90
- Passive mixing and collective Fock states: 349/349
- Loss and phase diffusion: 99/99
- Gaussian covariances and quadrature conventions: 169/169
- Photon-added centred Gaussian: 77/77
- Exact doubled positive-W: 140/140

## Other checks

- Full Aqua checks enabled, including undefined exports, compatibility, ambiguities and piracy.
- Documenter doctests passed; HTML documentation built with no warnings or errors. Local deployment disabled.
- All 23 executable state examples completed with an explicit seeded RNG.
- 66 warmed benchmark measurements completed for three batch sizes; all real occupation means within 2.45 estimated standard errors of their analytic targets.
- Independent finite-Fock/Kraus calculations and cutoff convergence are in the committed tests.
- [x,p]=i and [X,Y]=2i tested through X=√2x,Y=√2p, including vacuum covariances.
- `git diff --check` passed. GitHub-hosted CI has not been run/published from this local branch.
