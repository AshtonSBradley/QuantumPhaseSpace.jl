_modevalues(x::Number,M,name) = M == 1 ? [ComplexF64(x)] :
    throw(DimensionMismatch("$name requires $M entries."))
function _modevalues(x,M,name)
    length(x) == M || throw(DimensionMismatch("$name requires $M entries."))
    y = ComplexF64.(collect(x))
    all(isfinite,y) || throw(ArgumentError("$name must be finite."))
    return y
end

"""Q displacement `D(β) state D(β)†`. For pure states the overlap phase is retained."""
struct Displaced{S<:State} <: State
    state::S
    β::Vector{ComplexF64}
    function Displaced(s::S,β) where {S<:State}
        b = _modevalues(β,nmodes(s),"displacement")
        all(isfinite,b) || throw(ArgumentError("displacement must be finite."))
        new{S}(s,b)
    end
end
"""Rotate by `R(θ)=exp(iθ n)` per mode. A scalar rotates every mode equally."""
struct Rotated{S<:State} <: State
    state::S
    θ::Vector{Float64}
    function Rotated(s::S,θ) where {S<:State}
        angles = θ isa Real ? fill(Float64(θ),nmodes(s)) : Float64.(collect(θ))
        length(angles) == nmodes(s) || throw(DimensionMismatch("rotation requires one angle per mode."))
        all(isfinite,angles) || throw(ArgumentError("rotation angles must be finite."))
        new{S}(s,angles)
    end
end
nmodes(s::Union{Displaced,Rotated}) = nmodes(s.state)
overlapcapability(s::Union{Displaced,Rotated}) = overlapcapability(s.state)
function _qfill!(rng,Z,s::Displaced,d,limit)
    _qfill!(rng,Z,s.state,d,limit)
    Z .+= s.β
    return Z
end
function _qfill!(rng,Z,s::Rotated,d,limit)
    _qfill!(rng,Z,s.state,d,limit)
    Z .*= cis.(s.θ)
    return Z
end
_coordinate(z::Number) = [z]
_coordinate(z) = z
function logoverlap(s::Displaced,z)
    z = _coordinate(z)
    length(z) == nmodes(s) || throw(DimensionMismatch("overlap coordinate length."))
    w = z .- s.β
    return _logpoint(s.state,w) + sum((s.β.*conj.(z).-conj.(s.β).*z)./2)
end
function logoverlap(s::Rotated,z)
    z = _coordinate(z)
    length(z) == nmodes(s) || throw(DimensionMismatch("overlap coordinate length."))
    return _logpoint(s.state,z.*cis.(-s.θ))
end
_logpoint(s,z) = nmodes(s) == 1 ? logoverlap(s,only(z)) : logoverlap(s,z)

"""Independent tensor product of one or more states; component mode blocks retain their order."""
struct ProductState{T<:Tuple} <: State
    states::T
    function ProductState(states::T) where {T<:Tuple}
        !isempty(states) && all(s->s isa State,states) || throw(ArgumentError("product requires states."))
        new{T}(states)
    end
end
ProductState(states::State...) = ProductState(states)
ProductState(states::AbstractVector{<:State}) = ProductState(Tuple(states))
nmodes(s::ProductState) = sum(nmodes,s.states)
overlapcapability(s::ProductState) = all(x->overlapcapability(x) isa PureOverlap,s.states) ? PureOverlap() : SamplingOnly()
function _qfill!(rng,Z,s::ProductState,d,limit)
    offset = 0
    for state in s.states
        rows = offset+1:offset+nmodes(state)
        _qfill!(rng,view(Z,rows,:),state,SamplingDiagnostics(),limit)
        offset += nmodes(state)
    end
    d.proposals += size(Z,2)
    d.accepted += size(Z,2)
    return Z
end
function logoverlap(s::ProductState,z)
    length(z) == nmodes(s) || throw(DimensionMismatch("overlap coordinate length."))
    offset = 0
    result = 0.0im
    for state in s.states
        result += _logpoint(state,view(z,offset+1:offset+nmodes(state)))
        offset += nmodes(state)
    end
    return result
end

function _probabilities(p)
    q = Float64.(collect(p))
    !isempty(q) && all(x->isfinite(x) && x>=0,q) && maximum(q)>0 ||
        throw(ArgumentError("probabilities must be finite, nonnegative and not all zero."))
    q ./= maximum(q)
    q ./= sum(q)
    return q
end

"""Classical mixture with finite nonnegative weights (normalised on construction); never a coherent superposition."""
struct Mixture{T<:Tuple} <: State
    states::T
    probabilities::Vector{Float64}
    categorical::Categorical{Float64,Vector{Float64}}
    function Mixture(states,probabilities)
        states = Tuple(states)
        p = _probabilities(probabilities)
        length(states) == length(p) && all(s->s isa State,states) || throw(ArgumentError("one weight per state is required."))
        all(s->nmodes(s)==nmodes(first(states)),states) || throw(DimensionMismatch("mixture mode counts differ."))
        new{typeof(states)}(states,p,Categorical(p))
    end
end
nmodes(s::Mixture) = nmodes(first(s.states))
function _qfill!(rng,Z,s::Mixture,d,limit)
    components = rand(rng,s.categorical,size(Z,2))
    for j in eachindex(s.states)
        columns = findall(==(j),components)
        isempty(columns) && continue
        _qfill!(rng,view(Z,:,columns),s.states[j],SamplingDiagnostics(),limit)
    end
    d.proposals += size(Z,2)
    d.accepted += size(Z,2)
    return Z
end

"""Finite number-diagonal mixture. `p[n+1]` weights `|n⟩⟨n|`; no hidden cutoff is imposed."""
NumberMixture(p) = Mixture(Tuple(Fock(n) for n in 0:length(p)-1),p)

"""Normalised `(a†)^m ρthermal a^m`; nbar refers to the original thermal state."""
struct PhotonAddedThermal <: State
    m::Int
    nbar::Float64
    function PhotonAddedThermal(m::Integer,nbar::Real)
        m = _samplecount(m)
        isfinite(nbar) && nbar>=0 || throw(ArgumentError("nbar must be finite and nonnegative."))
        new(m,Float64(nbar))
    end
end
"""Normalised `a^m ρthermal (a†)^m`. For m>0, subtraction from vacuum is undefined."""
struct PhotonSubtractedThermal <: State
    m::Int
    nbar::Float64
    function PhotonSubtractedThermal(m::Integer,nbar::Real)
        m = _samplecount(m)
        isfinite(nbar) && nbar>=0 && (m==0 || nbar>0) || throw(ArgumentError("nbar must be finite, nonnegative, and positive when m>0."))
        new(m,Float64(nbar))
    end
end
function _qfill!(rng,Z,s::PhotonAddedThermal,d,limit)
    radial = Gamma(s.m+1,s.nbar+1)
    for i in eachindex(Z)
        Z[i] = sqrt(rand(rng,radial))*cis(2π*rand(rng))
    end
    d.proposals += length(Z)
    d.accepted += length(Z)
    return Z
end
function glauberP(rng::AbstractRNG,s::PhotonSubtractedThermal,N)
    N = _samplecount(N)
    a = zeros(ComplexF64,N)
    if s.nbar > 0
        radial = Gamma(s.m+1,s.nbar)
        for i in eachindex(a)
            a[i] = sqrt(rand(rng,radial))*cis(2π*rand(rng))
        end
    end
    return a,conj.(a)
end
positiveP(rng::AbstractRNG,s::PhotonSubtractedThermal,N) = glauberP(rng,s,N)
function wigner(rng::AbstractRNG,s::PhotonSubtractedThermal,N)
    a,_ = glauberP(rng,s,N)
    a .+= randnc(rng,length(a))/sqrt(2)
    return a,conj.(a)
end
function _qfill!(rng,Z,s::PhotonSubtractedThermal,d,limit)
    a,_ = glauberP(rng,s,size(Z,2))
    for i in eachindex(a)
        Z[1,i] = a[i] + randnc(rng)
    end
    d.proposals += length(a)
    d.accepted += length(a)
    return Z
end
