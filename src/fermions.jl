"""
    SlaterDeterminant(occupations)
    SlaterDeterminant(occupied_orbitals)

A number-conserving fermionic pure state. Supply a vector of integer occupations
(zero or one in the sampling basis), or an `M × N` matrix whose columns are
orthonormal occupied orbitals in an `M`-mode basis. Each mode is a spin-orbital.
"""
struct SlaterDeterminant
    occupations::Vector{Float64}
    orbitals::Union{Nothing,Matrix{ComplexF64}}

    function SlaterDeterminant(occupations::AbstractVector{<:Integer})
        !isempty(occupations) && all(n -> 0 ≤ n ≤ 1, occupations) ||
            throw(ArgumentError("supply a nonempty vector of zero or one occupations"))
        return new(Float64.(collect(occupations)), nothing)
    end

    function SlaterDeterminant(occupied_orbitals::AbstractMatrix{<:Number})
        M, N = size(occupied_orbitals)
        M > 0 && N ≤ M || throw(DimensionMismatch("expected M ≥ N occupied-orbital rows and M > 0"))
        C = Matrix{ComplexF64}(occupied_orbitals)
        all(isfinite, C) && isapprox(C' * C, I; atol=1e-10, rtol=1e-8) ||
            throw(ArgumentError("occupied orbitals must be finite and orthonormal"))
        occupations = vcat(ones(N), zeros(M-N))
        N in (0, M) && return new(occupations, nothing)
        # The implicit Q factor also gives an orthonormal basis for the empty orbitals.
        basis = lmul!(qr(C).Q, Matrix{ComplexF64}(I, M, M))
        return new(occupations, basis)
    end
end

struct CanonicalFermionOccupations
    particles::Int
    log_filled::Vector{Float64}
    log_partition::Matrix{Float64}
end

"""
    FermionThermal(energies, β; μ=nothing, particles=nothing)
    FermionThermal(h, β; μ=nothing, particles=nothing)

Free-fermion thermal state at finite inverse temperature `β ≥ 0`. Supply
one-body energies in the sampling basis or a Hermitian one-body Hamiltonian `h`.
By default the state is grand canonical with chemical potential zero. Set `μ`
for a different chemical potential, or set `particles` for the canonical
ensemble at that fixed total number.
"""
struct FermionThermal{T<:Union{Vector{Float64},CanonicalFermionOccupations}}
    occupations::T
    orbitals::Union{Nothing,Matrix{ComplexF64}}
end

function FermionThermal(energies::AbstractVector{<:Real}, β::Real; μ=nothing, particles=nothing)
    return _fermion_thermal(Float64.(collect(energies)), nothing, Float64(β), μ, particles)
end

function FermionThermal(h::AbstractMatrix{<:Number}, β::Real; μ=nothing, particles=nothing)
    size(h, 1) > 0 && ishermitian(h) && all(isfinite, h) ||
        throw(ArgumentError("the one-body Hamiltonian must be finite, nonempty and Hermitian"))
    H = Matrix{ComplexF64}(h)
    all(isfinite, H) || throw(ArgumentError("the one-body Hamiltonian must be finite in Float64"))
    (; values, vectors) = eigen(Hermitian(H))
    return _fermion_thermal(values, vectors, Float64(β), μ, particles)
end

function _fermion_thermal(energies, orbitals, β, μ, particles)
    M = length(energies)
    M > 0 && all(isfinite, energies) || throw(ArgumentError("supply finite one-body energies"))
    isfinite(β) && β ≥ 0 || throw(ArgumentError("β must be finite and nonnegative"))
    if isnothing(particles)
        chemical_potential = something(μ, 0.0)
        chemical_potential isa Real && isfinite(Float64(chemical_potential)) ||
            throw(ArgumentError("μ must be finite and real"))
        μ = Float64(chemical_potential)
        occupations = iszero(β) ? fill(0.5, M) : [inv(1 + exp(β*(ε-μ))) for ε in energies]
        return FermionThermal(occupations, orbitals)
    end
    isnothing(μ) || throw(ArgumentError("choose either a chemical potential or a fixed particle number"))
    particles isa Integer && 0 ≤ particles ≤ M ||
        throw(ArgumentError("particles must be an integer from zero to the number of modes"))
    particles in (0, M) && return FermionThermal(fill(Float64(particles == M), M), orbitals)
    return FermionThermal(_canonical_fermion_occupations(energies, β, Int(particles)), orbitals)
end

function _canonical_fermion_occupations(energies, β, particles)
    M = length(energies)
    reference = partialsort(energies, particles)
    log_filled = zeros(M)
    log_partition = fill(-Inf, M+1, particles+1)
    log_partition[1, 1] = 0
    # Shifting the spectrum by the Fermi energy leaves fixed-number weights
    # proportional. Per-mode normalization keeps every log weight nonpositive.
    for j in eachindex(energies)
        gap = iszero(β) ? 0.0 : β * (energies[j]-reference)
        log_filled[j] = -max(gap, 0.0)
        log_empty = min(gap, 0.0)
        log_partition[j+1, 1] = log_partition[j, 1] + log_empty
        for k in 1:min(j, particles)
            empty = log_partition[j, k+1] + log_empty
            filled = log_partition[j, k] + log_filled[j]
            hi, lo = max(empty, filled), min(empty, filled)
            log_partition[j+1, k+1] = isfinite(hi) ? hi + log1p(exp(lo-hi)) : hi
        end
    end
    return CanonicalFermionOccupations(particles, log_filled, log_partition)
end

function _draw_fermion_occupations!(rng, occupations, state::CanonicalFermionOccupations)
    (; particles, log_filled, log_partition) = state
    fill!(occupations, 0)
    remaining = particles
    for j in length(occupations):-1:1
        remaining == 0 && break
        log_p = log_filled[j] + log_partition[j, remaining] - log_partition[j+1, remaining+1]
        if log(rand(rng)) < log_p
            occupations[j] = 1
            remaining -= 1
        end
    end
    return occupations
end

nmodes(s::SlaterDeterminant) = length(s.occupations)
nmodes(s::FermionThermal{Vector{Float64}}) = length(s.occupations)
nmodes(s::FermionThermal{CanonicalFermionOccupations}) = length(s.occupations.log_filled)

"""
    wigner([rng], state::Union{SlaterDeterminant,FermionThermal}, N; pairing=false)

Moment-matched initial ensemble for fermionic truncated Wigner. Return a
complex `M × M × N` array `ρ`, where `ρ[i,j,k]` is a sample of
`cᵢ†cⱼ - δᵢⱼ/2`. With `pairing=true`, return `(; ρ, τ)` with the antisymmetric
pair variables `τ[i,j,k]` for `cᵢcⱼ`. The Gaussian draw matches the means and
symmetrized covariances of the bilinears. For a fixed-number thermal state,
first draw a Slater determinant from the canonical ensemble and apply the
moment matching conditional on that determinant.
"""
function wigner(rng::AbstractRNG, state::Union{SlaterDeterminant,FermionThermal}, N; pairing::Bool=false)
    M, N = nmodes(state), _samplecount(N)
    ρ = Array{ComplexF64}(undef, M, M, N)
    τ = pairing ? similar(ρ) : nothing
    N > 0 && _fermion_wigner!(rng, ρ, τ, state.occupations, state.orbitals)
    return pairing ? (; ρ, τ) : ρ
end

wigner(state::Union{SlaterDeterminant,FermionThermal}, N; kwargs...) =
    wigner(Random.default_rng(), state, N; kwargs...)

function _fermion_density!(rng, ρ, occupations)
    for j in eachindex(occupations)
        fj = occupations[j]
        vacancy = 1-fj
        ρ[j,j] = fj - 0.5 + (iszero(fj*vacancy) ? 0.0 : sqrt(fj*vacancy)*randn(rng))
        for i in 1:j-1
            fi = occupations[i]
            variance = (fi*vacancy + fj*(1-fi))/2
            value = iszero(variance) ? 0.0im : sqrt(variance)*randnc(rng)
            ρ[i,j] = value
            ρ[j,i] = conj(value)
        end
    end
    return ρ
end

function _fermion_pairs!(rng, τ, occupations)
    for j in eachindex(occupations)
        fj = occupations[j]
        τ[j,j] = 0
        for i in 1:j-1
            fi = occupations[i]
            variance = (fi*fj + (1-fi)*(1-fj))/2
            value = iszero(variance) ? 0.0im : sqrt(variance)*randnc(rng)
            τ[i,j] = value
            τ[j,i] = -value
        end
    end
    return τ
end

function _fermion_wigner!(rng, ρ, τ, occupation_state, orbitals)
    M = size(ρ, 1)
    occupations = occupation_state isa CanonicalFermionOccupations ? zeros(M) : occupation_state
    rotate = !isnothing(orbitals)
    natural_density = rotate ? Matrix{ComplexF64}(undef, M, M) : nothing
    natural_pairs = rotate && !isnothing(τ) ? similar(natural_density) : nothing
    work = rotate ? similar(natural_density) : nothing
    conjugate_orbitals = rotate ? conj.(orbitals) : nothing
    for k in axes(ρ, 3)
        if occupation_state isa CanonicalFermionOccupations
            _draw_fermion_occupations!(rng, occupations, occupation_state)
        end
        density = view(ρ, :, :, k)
        _fermion_density!(rng, rotate ? natural_density : density, occupations)
        if rotate
            mul!(work, conjugate_orbitals, natural_density)
            mul!(density, work, transpose(orbitals))
            hermitianpart!(density)
        end
        if !isnothing(τ)
            pairs = view(τ, :, :, k)
            _fermion_pairs!(rng, rotate ? natural_pairs : pairs, occupations)
            if rotate
                mul!(work, orbitals, natural_pairs)
                mul!(pairs, work, transpose(orbitals))
                for j in 1:M
                    pairs[j,j] = 0
                    for i in 1:j-1
                        value = (pairs[i,j]-pairs[j,i])/2
                        pairs[i,j] = value
                        pairs[j,i] = -value
                    end
                end
            end
        end
    end
    return ρ
end
