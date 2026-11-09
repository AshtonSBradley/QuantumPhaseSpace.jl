# PhaseSpaceTools

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://AshtonSBradley.github.io/PhaseSpaceTools.jl/dev)
[![Build Status](https://github.com/AshtonSBradley/PhaseSpaceTools.jl/workflows/CI/badge.svg)](https://github.com/AshtonSBradley/PhaseSpaceTools.jl/actions)
[![Aqua](https://github.com/AshtonSBradley/PhaseSpaceTools.jl/actions/workflows/Aqua.yml/badge.svg)](https://github.com/AshtonSBradley/PhaseSpaceTools.jl/actions/workflows/Aqua.yml)
[![Coverage](https://codecov.io/gh/AshtonSBradley/PhaseSpaceTools.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/AshtonSBradley/PhaseSpaceTools.jl)
[![DOI](https://zenodo.org/badge/115932136.svg)](https://zenodo.org/badge/latestdoi/115932136)

Sample quantum initial states commonly encountered in quantum phase space simulations of Bose fields, including those encountered in quantum optics and Bose-Einstein condensates. 

Wigner and positive-P distributions are available, being the most useful for dynamical simulations.

Available distributions are `glauberP`, `positiveP`, `wigner`, `positiveW`, and `husimiQ`.

<img src="examples/fock.png" alt="Fock state Wigner function" width="300"/>

## Install

```julia
julia> ]add PhaseSpaceTools
```

## Quick Start

```julia
julia> using PhaseSpaceTools

julia> α = 1.0 + 2.0im
julia> state = Coherent(α)
julia> samples, samples_conj = positiveP(state, 4)

julia> mean(samples)
1.0 + 2.0im

julia> real(mean(samples_conj .* samples))
5.0
```

For coherent states, the `+P` representation is deterministic: `samples` and `samples_conj` are fixed at `α` and `conj(α)`.

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

The [support table](docs/src/support.md) distinguishes exact laws, explicit
approximations and unsupported methods. Positive-W is a doubled representation,
with generally correlated nonconjugate coordinates, and is not ordinary W.
Gaussian covariances use `[x,p]=i` and vacuum `I/2`; divide an `[X,Y]=2i`
covariance by two before using `Gaussian`.

Development evidence and outstanding limitations are recorded in
[`dev/PROGRESS.md`](dev/PROGRESS.md).

# External links
___Numerical representation of quantum states in the positive-P and Wigner representations,___ \
M K Olsen, A S Bradley, \
[Optics Communications 282, 3924 (2009)](http://dx.doi.org/10.1016/j.optcom.2009.06.033)

[Scipost Commentary with full erratum](https://scipost.org/commentaries/10.1016/j.optcom.2009.06.033/)
