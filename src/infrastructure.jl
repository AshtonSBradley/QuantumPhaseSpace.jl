_samplecount(N) = N isa Integer && 0 <= N <= typemax(Int) ? Int(N) :
    throw(ArgumentError("sample/attempt count must be a nonnegative integer representable as Int."))
_attemptlimit(N) = max(1000,Base.checked_mul(1000,max(_samplecount(N),1)))

"""Number of modes. Matrix sampling uses modes × samples."""
nmodes(::State) = 1
nmodes(::SqueezedTwoMode) = 2

"""Proposal and accepted counts for the outer rejection sampler (zero for an empty batch)."""
mutable struct SamplingDiagnostics
    proposals::Int
    accepted::Int
end
SamplingDiagnostics() = SamplingDiagnostics(0,0)
"""Observed acceptance; NaN if no proposals were required."""
acceptance_rate(d::SamplingDiagnostics) = d.proposals == 0 ? NaN : d.accepted/d.proposals

"""Attempt-limit exhaustion. `diagnostics` records progress; no partial batch is returned."""
struct SamplingError <: Exception
    limit::Int
    diagnostics::SamplingDiagnostics
end
Base.showerror(io::IO,e::SamplingError) = print(io,"sampling exhausted ",e.limit,
    " proposals; accepted ",e.diagnostics.accepted," (rate ",acceptance_rate(e.diagnostics),
    "). Increase max_attempts or choose a more efficient exact proposal.")
function _proposal!(d,limit)
    d.proposals < limit || throw(SamplingError(limit,deepcopy(d)))
    d.proposals += 1
end

abstract type OverlapCapability end
"""A normalised pure state with evaluable complex coherent-state overlaps."""
struct PureOverlap <: OverlapCapability end
"""Sampling capability alone does not supply a pure-state overlap."""
struct SamplingOnly <: OverlapCapability end
"""Return `PureOverlap()` or `SamplingOnly()`; mixed states cannot be superposition components."""
overlapcapability(::State) = SamplingOnly()
overlapcapability(::Union{Coherent,Fock,Squeezed}) = PureOverlap()

"""
    logoverlap(state, z)

Complex logarithm of ⟨z|state⟩, including its phase. The real part is the log
magnitude. Zeros have real part -Inf; sampling engines cancel common scales.
"""
logoverlap(s::Coherent,z::Number) = -abs2(z-s.β)/2 + im*imag(conj(z)*s.β)
function logoverlap(s::Fock,z::Number)
    s.n == 0 && return complex(-abs2(z)/2)
    iszero(z) && return complex(-Inf)
    return complex(-abs2(z)/2 + s.n*log(abs(z)) - loggamma(s.n+1)/2,
                   -s.n*angle(z))
end
function logoverlap(s::Squeezed,z::Number)
    (; β,ϵ) = s
    r = abs(ϵ)
    w = z-β
    t = -cis(angle(ϵ))*tanh(r)
    return -abs2(w)/2 + t*conj(w)^2/2 - log(cosh(r))/2 +
        (β*conj(z)-conj(β)*z)/2
end

_legacy(s,A,B) = nmodes(s) == 1 ? (vec(A),vec(B)) : (A,B)
_legacy(::SqueezedTwoMode,A,B) = (vec(A[1,:]),vec(B[1,:]),vec(A[2,:]),vec(B[2,:]))

"""
    sampleQ([rng], state, N; diagnostics=false, max_attempts=max(1000,1000N))

Return one complex matrix of Q coordinates, with modes in rows and samples
in columns, for every state (including legacy two-mode states). With
`diagnostics=true`, return `(matrix, SamplingDiagnostics)`. Rejection limits
apply to outer proposals; nested component samplers enforce their own limits.
"""
function sampleQ(rng::AbstractRNG,s::State,N;diagnostics=false,max_attempts=_attemptlimit(N))
    Z = Matrix{ComplexF64}(undef,nmodes(s),_samplecount(N))
    d = SamplingDiagnostics()
    _qfill!(rng,Z,s,d,_samplecount(max_attempts))
    return diagnostics ? (Z,d) : Z
end
sampleQ(s::State,N;kwargs...) = sampleQ(Random.default_rng(),s,N;kwargs...)

"""
    sampleQ!([rng], Z, state; max_attempts=...)

Fill a caller-owned complex modes × samples matrix. Returns `Z`. On rejection
failure the destination may contain a partial batch and a `SamplingError` is
raised; only successful returns represent complete batches.
"""
function sampleQ!(rng::AbstractRNG,Z::AbstractMatrix{<:Complex},s::State;
                  max_attempts=_attemptlimit(size(Z,2)))
    size(Z,1) == nmodes(s) || throw(DimensionMismatch("expected $(nmodes(s)) mode rows."))
    _qfill!(rng,Z,s,SamplingDiagnostics(),_samplecount(max_attempts))
    return Z
end
sampleQ!(Z,s::State;kwargs...) = sampleQ!(Random.default_rng(),Z,s;kwargs...)

const LegacyQState = Union{Coherent,Fock,Thermal,Squeezed,Crescent,Bogoliubov,SqueezedTwoMode}
function _qfill!(rng,Z,s::LegacyQState,d,limit)
    result = husimiQ(rng,s,size(Z,2))
    for m in 1:nmodes(s)
        Z[m,:] .= result[2m-1]
    end
    d.proposals += size(Z,2)
    d.accepted += size(Z,2)
    return Z
end
function husimiQ(rng::AbstractRNG,s::State,N;kwargs...)
    Z = sampleQ(rng,s,N;kwargs...)
    Z isa Tuple && throw(ArgumentError("use sampleQ(...;diagnostics=true) for diagnostics."))
    return _legacy(s,Z,conj.(Z))
end
