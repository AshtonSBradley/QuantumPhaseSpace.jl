"""
    Passive(state, U)

Passive mixing `a_out=U*a_in`, with a square unitary matrix. Include vacuum
input modes explicitly. Q coordinates transform by U. A structured operator
can extend `_qfill!` using the same `mul!` contract; the public constructor
currently validates and caches an explicit matrix.
"""
struct Passive{S<:State} <: State
    state::S
    U::Matrix{ComplexF64}
    function Passive(s::S,U::AbstractMatrix) where {S<:State}
        size(U)==(nmodes(s),nmodes(s)) || throw(DimensionMismatch("U must match input mode count."))
        U = Matrix{ComplexF64}(U)
        all(isfinite,U) && isapprox(U'*U,I;atol=1e-10,rtol=1e-10) || throw(ArgumentError("U must be finite and unitary."))
        new{S}(s,U)
    end
end
nmodes(s::Passive) = nmodes(s.state)
overlapcapability(s::Passive) = overlapcapability(s.state)
function logoverlap(s::Passive,z)
    length(z)==nmodes(s) || throw(DimensionMismatch("overlap coordinate length."))
    return _logpoint(s.state,s.U'*z)
end
function _qfill!(rng,Z,s::Passive,d,limit)
    input = similar(Z)
    _qfill!(rng,input,s.state,d,limit)
    mul!(Z,s.U,input)
    return Z
end

"""
    FixedNumberCondensate(N, u)

State `(Σ u[j]a[j]†)^N|0⟩/√N!`, with `u'u=1`. The exact O(M) Q sampler
uses a Fock collective mode and independent vacuum in its orthogonal complement.
Sample norms fluctuate and are never renormalised. `positiveP` defaults to the
specialised occupied-mode construction; `method=:generic` selects Q conversion.
"""
struct FixedNumberCondensate <: State
    N::Int
    u::Vector{ComplexF64}
    function FixedNumberCondensate(N::Integer,u::AbstractVector)
        N = _samplecount(N)
        u = ComplexF64.(collect(u))
        !isempty(u) && all(isfinite,u) && isapprox(sum(abs2,u),1;atol=1e-12,rtol=1e-12) ||
            throw(ArgumentError("u must be finite and normalised (u†u=1)."))
        new(N,u)
    end
end
nmodes(s::FixedNumberCondensate) = length(s.u)
overlapcapability(::FixedNumberCondensate) = PureOverlap()
function logoverlap(s::FixedNumberCondensate,z)
    z = _coordinate(z)
    length(z)==nmodes(s) || throw(DimensionMismatch("overlap coordinate length."))
    w = dot(s.u,z)
    return logoverlap(Fock(s.N),w)-(sum(abs2,z)-abs2(w))/2
end
function _qfill!(rng,Z,s::FixedNumberCondensate,d,limit)
    (; N,u) = s
    radial = Gamma(N+1,1)
    for col in eachcol(Z)
        randn!(rng,col)
        projection = dot(u,col)
        occupied = sqrt(rand(rng,radial))*cis(2π*rand(rng))
        col .+= u.*(occupied-projection)
    end
    d.proposals += size(Z,2)
    d.accepted += size(Z,2)
    return Z
end
function positiveP(rng::AbstractRNG,s::FixedNumberCondensate,N;method=:specialised)
    N = _samplecount(N)
    if method == :generic
        μ = sampleQ(rng,s,N)
        γ = randnc(rng,size(μ)...)
        return _legacy(s,μ.+γ,conj.(μ.-γ))
    elseif method != :specialised
        throw(ArgumentError("method must be :specialised or :generic."))
    end
    # Deterministic positive-P vacuum in every orthogonal collective mode.
    if s.N == 0
        Z = zeros(ComplexF64,nmodes(s),N)
        return _legacy(s,Z,copy(Z))
    end
    a,ap = positiveP(rng,Fock(s.N),N)
    return _legacy(s,s.u*transpose(a),conj.(s.u)*transpose(ap))
end
