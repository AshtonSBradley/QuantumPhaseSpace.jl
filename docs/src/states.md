# Quantum states
Quantum states available in `PhaseSpaceTools.jl` are

## Coherent
```@docs
Coherent
```
## Fock
```@docs
Fock
```
## Crescent
```@docs
Crescent
```
## Squeezed
```@docs
Squeezed
```
## SqueezedTwoMode
```@docs
SqueezedTwoMode
```
## Thermal
```@docs
Thermal
```
## Bogoliubov
```@docs
Bogoliubov
```

## Cats and coherent superpositions

[`EvenCat`](@ref), [`OddCat`](@ref) and [`Cat`](@ref) combine coherent components at
±α with a chosen relative phase. A [`Mixture`](@ref) of those components is a
different density matrix. Follow the [Q sampling tutorial](@ref cat-q-tutorial)
for figures, runnable code and moment checks; the
[finite-superposition machinery](exact-sampling.md) also covers sparse Fock
states, NOON states and photon-added coherent states. The
[Fock-basis examples](@ref fock-basis-examples) show how to check a
number-state coherence and sample a Kerr-evolved state.

## Exact Q compositions and thermal operations

Displacement and rotation act on Q coordinates as `z+β` and `exp(iθ)z`.
`Mixture` samples a density-matrix mixture, whereas `ProductState` concatenates
independent mode blocks. Neither operation introduces coherent interference.
Finite number mixtures are exact for the supplied probabilities; if these
were obtained by truncation, the caller must report that cutoff and tail mass.

Photon-added thermal Q intensities have Gamma shape `m+1` and scale `nbar+1`.
Photon-subtracted thermal P intensities have shape `m+1` and scale `nbar`;
adding independent circular noise gives W (variance 1/2) or Q (variance 1).
These Gamma parameters always denote shape and scale.

```@docs
Displaced
Rotated
Mixture
ProductState
NumberMixture
PhotonAddedThermal
PhotonSubtractedThermal
```
