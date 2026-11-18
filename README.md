# QuantumPhaseSpace

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ashtonsbradley.github.io/QuantumPhaseSpace.jl/)
[![Build Status](https://github.com/AshtonSBradley/QuantumPhaseSpace.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/AshtonSBradley/QuantumPhaseSpace.jl/actions/workflows/CI.yml)
[![Aqua](https://github.com/AshtonSBradley/QuantumPhaseSpace.jl/actions/workflows/Aqua.yml/badge.svg?branch=master)](https://github.com/AshtonSBradley/QuantumPhaseSpace.jl/actions/workflows/Aqua.yml)
[![Coverage](https://codecov.io/gh/AshtonSBradley/QuantumPhaseSpace.jl/branch/master/graph/badge.svg)](https://app.codecov.io/gh/AshtonSBradley/QuantumPhaseSpace.jl)
[![PhaseSpaceTools DOI](https://img.shields.io/badge/legacy%20DOI-10.5281%2Fzenodo.3758700-blue)](https://doi.org/10.5281/zenodo.3758700)

Sample initial states for bosonic and fermionic quantum phase-space simulations.
The bosonic samplers cover quantum optics and Bose fields. The fermionic samplers
provide initial bilinear ensembles for truncated Wigner calculations.

QuantumPhaseSpace.jl succeeds [PhaseSpaceTools.jl](https://github.com/AshtonSBradley/PhaseSpaceTools.jl).
See the [documentation](docs/src/index.md) for usage and representation conventions.

Exact Husimi-Q and positive-P sampling includes cats, finite quantum superpositions, photon-added states and multimode initial conditions. Ordinary Wigner sampling is available for the states listed in the support table.

Available distributions are `glauberP`, `positiveP`, `wigner`, `positiveW`, and `husimiQ`.

[![Even-cat, odd-cat and classical-mixture Husimi-Q samples, with coherent centres at plus and minus alpha.](docs/src/assets/cat-q.png)](docs/src/cat-q.md)

Even and odd cats preserve different interference between the same coherent
components at ±α. [Explore the cat-state Q tutorial](docs/src/cat-q.md).

## Install

```julia
using Pkg
Pkg.add(url="https://github.com/AshtonSBradley/QuantumPhaseSpace.jl")
```

## Quick Start

```julia
using QuantumPhaseSpace, Random

rng = MersenneTwister(42)
α = 1.2 + 0im
state = EvenCat(α)                  # |α⟩ + |−α⟩, normalised
z, zbar = husimiQ(rng, state, 100_000)
```

`z` contains Q phase-plane samples and `zbar == conj.(z)`. Use `OddCat(α)` for
the relative minus sign. The classical mixture is
`Mixture((Coherent(α),Coherent(-α)),[0.5,0.5])`. The
[tutorial](docs/src/cat-q.md) compares all three distributions and explains
`sampleQ(...; diagnostics=true)`, ordering corrections and preparation loss.

## Core states

- `Bogoliubov`
- `Coherent`
- `Crescent`
- `Fock`
- `Squeezed`
- `SqueezedTwoMode`
- `Thermal`

Exact Q extensions include cats and finite coherent/Fock superpositions, NOON
states, photon-added and photon-subtracted families, passive multimode states,
fixed-number condensates, Gaussian covariance states, and loss/phase-diffusion
channels. All inherit exact positive-P sampling. See the
[state × representation table](docs/src/support.md) for the full API and its
exact/approximate/unsupported distinctions.

### Coherent state
A coherent state |α⟩ is sampled as
```julia
α = 1.0+im*2.0 # coherent amplitude
s = Coherent(α) # define state |α⟩
N = 1000 # number of samples
a,a⁺ = positiveP(s,N)
```
This is a special case where the two phase space variables `a` and `a⁺` are complex conjugate, and non-stochastic in the `+P` representation.

### Fock state

<img src="examples/fock.png" alt="Fock state Wigner function" width="300"/>

An approximate Fock state sampler in the Wigner representation:
```julia
n = 100
s = Fock(n) # define number state |n⟩
N = 1000 # number of samples
a,a⁺ = wigner(s,N)
```
Provides an approximate sampling of `W` that reproduces operator averages for large `n`. For every `n>0` this is explicitly an approximation, useful when `n ≫ 1`. `Fock(0)` uses the exact vacuum Gaussian.

## Examples

See [`examples/exact_sampling.jl`](examples/exact_sampling.jl) for seeded examples
of cats, photon addition, NOON states, fixed-number condensates and preparation
loss. All samplers accept an explicit RNG as their first argument. `sampleQ`
always returns a modes × samples matrix; historical single-mode vectors and the
legacy two-mode four-vector interfaces remain available.

- [Kerr evolution](examples/kerr_evolution.jl) samples a coherent input at
  collapse and revival times using finite Fock superpositions.
- [Canonical thermal states](examples/canonical_thermal.jl) sample Gibbs
  mixtures for Kerr, harmonic and displaced harmonic oscillators.
- The [worked Fock-basis examples](docs/src/fock-basis.md) explain the
  estimators and checks used in these examples.
- [Fermionic truncated Wigner](docs/src/fermions.md) initializes bilinear
  trajectories for Slater determinants and free thermal fermions.

The [support table](docs/src/support.md) distinguishes exact laws, explicit
approximations and unsupported methods. Positive-W is a doubled representation,
with generally correlated nonconjugate coordinates, and is not ordinary W.
Gaussian covariances use `[x,p]=i` and vacuum `I/2`; divide an `[X,Y]=2i`
covariance by two before using `Gaussian`.

Development evidence and outstanding limitations are recorded in
[`dev/PROGRESS.md`](dev/PROGRESS.md).

## References

- M. K. Olsen and A. S. Bradley, [Numerical representation of quantum states in the positive-P and Wigner representations](https://doi.org/10.1016/j.optcom.2009.06.033), *Optics Communications* **282**, 3924–3929 (2009). See the [2016 thermal-state erratum](https://doi.org/10.1016/j.optcom.2016.02.068) and the [SciPost commentary](https://scipost.org/commentaries/10.1016/j.optcom.2009.06.033/) for corrections to the original prescriptions.
- P. D. Drummond and C. W. Gardiner, [Generalised P-representations in quantum optics](https://doi.org/10.1088/0305-4470/13/7/018), *J. Phys. A* **13**, 2353–2368 (1980). This gives the canonical positive-P construction used with the package's Husimi-Q samplers.
- P. D. Drummond and B. Opanchuk, [Initial states for quantum field simulations in phase space](https://doi.org/10.1103/PhysRevResearch.2.033304), *Phys. Rev. Research* **2**, 033304 (2020), and P. D. Drummond et al., [Simulating complex networks in phase space: Gaussian boson sampling](https://doi.org/10.1103/PhysRevA.105.012427), *Phys. Rev. A* **105**, 012427 (2022). These provide the context for direct Gaussian positive-P sampling.
- L. I. Plimak et al., [Beyond the Fokker–Planck equation: Stochastic simulation of complete Wigner representation for the optical parametric oscillator](https://doi.org/10.1209/epl/i2001-00529-8), *Europhys. Lett.* **56**, 372–378 (2001). This develops the doubled positive-W representation.
- S. M. Davidson, D. Sels and A. Polkovnikov, [Semiclassical approach to dynamics of interacting fermions](https://doi.org/10.1016/j.aop.2017.07.003), *Ann. Phys.* **384**, 128–141 (2017). This gives the fermionic truncated-Wigner framework used by the bilinear initializers.
