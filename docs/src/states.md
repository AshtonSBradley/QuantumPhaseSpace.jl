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
