# Overview
This package provides sampling methods for commonly used quantum states in various quantum phase-space representations. These include Glauber–Sudarshan P, Husimi Q, ordinary Wigner, and the doubled positive-P and positive-W representations.

It also provides helpers for sampling complex and real noises used in phase-space simulations, including SDE workflows with [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl).

## Rejection sampling in action

[Exact Q Sampling](exact-sampling.md) introduces a rejection method for coherent
superpositions. The following [cat-state example](@ref cat-q-tutorial) shows how
acceptance creates the interference between components centred at +α and −α.
It compares even and odd cats with a classical mixture, then checks moments
and acceptance using seeded samples.

![Even-cat, odd-cat and classical-mixture Q samples at alpha=1.2, with exact-density contours.](assets/cat-q.png)

The [Fock-basis examples](@ref fock-basis-examples) use the same sampling
machinery for number-state coherences, coherent input evolved under a Kerr
Hamiltonian, and canonical thermal ensembles. The [support table](support.md)
lists the available samplers for each state.

# Installation
In the julia REPL

```julia-repl
julia> ]add PhaseSpaceTools

julia> ]test PhaseSpaceTools
```
# References
Cite [PhaseSpaceTools.jl](https://doi.org/10.5281/zenodo.3758700) for the
software. For the original sampling prescriptions, see

M. K. Olsen, A. S. Bradley, Numerical representation of quantum states in the positive-P and Wigner representations, [Optics Communications __282__, 3924-3929 (2009).](https://dx.doi.org/10.1016/j.optcom.2009.06.033)

See the [published thermal-state erratum](https://doi.org/10.1016/j.optcom.2016.02.068)
and the [SciPost commentary](https://scipost.org/commentaries/10.1016/j.optcom.2009.06.033/)
for corrections to the sampling prescriptions.
