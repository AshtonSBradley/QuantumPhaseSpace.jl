function _realmodevalues(x,M,name;lower=0.0,upper=Inf)
    y = x isa Real ? fill(Float64(x),M) : Float64.(collect(x))
    length(y)==M || throw(DimensionMismatch("$name requires $M values."))
    all(t->isfinite(t) && lower<=t<=upper,y) || throw(ArgumentError("$name must be finite and in [$lower,$upper]."))
    return y
end
"""
    Loss(state, η)

Exact pure-loss quantum channel on Q samples. Transmissivity is a scalar or
one value per mode, in [0,1]. Independent vacuum environment noise is included;
loss of a superposition includes decoherence and is not a pure shrunken cat.
"""
struct Loss{S<:State} <: State
    state::S
    η::Vector{Float64}
    function Loss(s::S,η) where {S<:State}
        new{S}(s,_realmodevalues(η,nmodes(s),"transmissivity";upper=1.0))
    end
end
"""
    PhaseDiffused(state, σ; common=true)

Gaussian phase diffusion. With `common=true`, a scalar standard deviation
rotates all modes by the same random angle. With `common=false`, a scalar or
per-mode vector specifies independent mode phases. σ must be nonnegative.
"""
struct PhaseDiffused{S<:State} <: State
    state::S
    σ::Vector{Float64}
    common::Bool
    function PhaseDiffused(s::S,σ;common::Bool=true) where {S<:State}
        common && !(σ isa Real) && throw(ArgumentError("common phase diffusion requires a scalar σ."))
        new{S}(s,_realmodevalues(σ,nmodes(s),"phase standard deviation"),common)
    end
end
nmodes(s::Union{Loss,PhaseDiffused}) = nmodes(s.state)
function _qfill!(rng,Z,s::Loss,d,limit)
    _qfill!(rng,Z,s.state,d,limit)
    for col in eachcol(Z), m in eachindex(s.η)
        η = s.η[m]
        η==1 && continue
        col[m] = sqrt(η)*col[m]+sqrt(1-η)*randnc(rng)
    end
    return Z
end
function _qfill!(rng,Z,s::PhaseDiffused,d,limit)
    _qfill!(rng,Z,s.state,d,limit)
    all(iszero,s.σ) && return Z
    for col in eachcol(Z)
        if s.common
            col .*= cis(s.σ[1]*randn(rng))
        else
            for m in eachindex(s.σ)
                col[m] *= cis(s.σ[m]*randn(rng))
            end
        end
    end
    return Z
end
