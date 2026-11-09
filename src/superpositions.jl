_primitives(s::State) = throw(ArgumentError("analytic pairwise Gram preparation is not implemented for $(typeof(s)); use the documented coherent/Fock/squeezed components or expand their sums componentwise."))

# Normalised Bargmann primitives: displaced Fock (κ=0), or pure Gaussian (n=0).
_primitive(s::Coherent) = (β=s.β,κ=0.0im,n=0,phase=0.0)
_primitive(s::Fock) = (β=0.0im,κ=0.0im,n=s.n,phase=0.0)
_primitive(s::Squeezed) = (β=s.β,κ=-cis(angle(s.ϵ))*tanh(abs(s.ϵ)),n=0,phase=0.0)
_primitives(s::Union{Coherent,Fock,Squeezed}) = [_primitive(s)]
_primitives(s::ProductState) = reduce(vcat,(_primitives(x) for x in s.states))
function _primitives(s::Displaced)
    return [merge(p,(β=p.β+b,phase=p.phase+imag(b*conj(p.β))))
            for (p,b) in zip(_primitives(s.state),s.β)]
end
function _primitives(s::Rotated)
    return [merge(p,(β=p.β*cis(θ),κ=p.κ*cis(2θ),phase=p.phase+p.n*θ))
            for (p,θ) in zip(_primitives(s.state),s.θ)]
end

# Independent Bargmann Gaussian integration, with normalised polynomial moments.
# Preparation uses BigFloat to avoid cancellation in displaced number overlaps.
function _primitive_inner(p,q)
    b,c = Complex{BigFloat}(p.β),Complex{BigFloat}(q.β)
    k,l = Complex{BigFloat}(p.κ),Complex{BigFloat}(q.κ)
    if iszero(k) && iszero(l) && b == c
        return p.n == q.n ? cis(BigFloat(q.phase-p.phase)) : zero(b)
    end
    h,j = b-k*conj(b),c-l*conj(c)
    C = exp(-abs2(b)/2+k*conj(b)^2/2)*(1-abs2(k))^(BigFloat(1)/4)
    D = exp(-abs2(c)/2+l*conj(c)^2/2)*(1-abs2(l))^(BigFloat(1)/4)
    den = 1-conj(k)*l
    base = conj(C)*D/sqrt(den)*exp((conj(h)*j+(conj(k)*j^2+l*conj(h)^2)/2)/den)
    ma = (j+l*conj(h))/den-b
    mb = (conj(h)+conj(k)*j)/den-conj(c)
    ca,cb,cab = l/den,conj(k)/den,1/den
    M = zeros(Complex{BigFloat},p.n+1,q.n+1)
    M[1,1] = 1
    for t in 1:q.n
        M[1,t+1] = mb/sqrt(BigFloat(t))*M[1,t]
        t>1 && (M[1,t+1] += cb*sqrt(BigFloat(t-1)/t)*M[1,t-1])
    end
    for r in 1:p.n, t in 0:q.n
        M[r+1,t+1] = ma/sqrt(BigFloat(r))*M[r,t+1]
        r>1 && (M[r+1,t+1] += ca*sqrt(BigFloat(r-1)/r)*M[r-1,t+1])
        t>0 && (M[r+1,t+1] += cab*sqrt(BigFloat(t)/r)*M[r,t])
    end
    return base*M[end,end]*cis(BigFloat(q.phase-p.phase))
end
function _stateinner(s::State,t::State)
    nmodes(s) == nmodes(t) || throw(DimensionMismatch("overlap mode counts differ."))
    return prod(_primitive_inner(p,q) for (p,q) in zip(_primitives(s),_primitives(t)))
end

"""
    Superposition(states, coefficients)

Normalised coherent sum of normalised pure states. Caches the Gram matrix and
proposal categorical law proportional to `abs.(coefficients)`. Supports coherent,
Fock, pure squeezed states, their products/displacements/rotations, and nested
superpositions. No Hilbert-space truncation is used. `expected_acceptance` is
Z/(Σ|c|)². Preparation rejects zero or numerically unresolved states.
"""
struct Superposition{T<:Tuple} <: State
    states::T
    coefficients::Vector{ComplexF64}
    gram::Matrix{ComplexF64}
    normalisation::Float64
    categorical::Categorical{Float64,Vector{Float64}}
    expected_acceptance::Float64
    function Superposition(states,coefficients)
        states = Tuple(states)
        c = ComplexF64.(collect(coefficients))
        !isempty(c) && length(c)==length(states) && all(isfinite,c) && maximum(abs,c)>0 ||
            throw(ArgumentError("provide one finite coefficient per state, with a nonzero coefficient."))
        all(s->s isa State && overlapcapability(s) isa PureOverlap,states) ||
            throw(ArgumentError("superposition requires pure states with complex overlaps, not mixed-state samplers."))
        all(s->nmodes(s)==nmodes(first(states)),states) || throw(DimensionMismatch("component mode counts differ."))
        keep = findall(!iszero,c)
        states,c = states[keep],c[keep]
        c ./= maximum(abs,c)
        c ./= sum(abs,c)
        precision = max(256,128+4maximum(p.n for s in states for p in _preparation_primitives(s)))
        G,Z = setprecision(BigFloat,precision) do
            g = [_stateinner(s,t) for s in states,t in states]
            bc = Complex{BigFloat}.(c)
            z = dot(bc,g*bc)
            isfinite(z) && abs(imag(z)) < big"1e-30" && real(z)>0 ||
                throw(ArgumentError("superposition is zero or its Gram normalisation is unresolved."))
            Matrix{ComplexF64}(g),Float64(real(z))
        end
        isfinite(Z) && Z>0 || throw(ArgumentError("normalisation is below Float64 range."))
        minimum(eigvals(Hermitian(G))) >= -1e-10 || throw(ArgumentError("component Gram matrix is not positive semidefinite."))
        p = abs.(c)
        new{typeof(states)}(states,c,G,Z,Categorical(p),Z)
    end
end
_preparation_primitives(s::State) = _primitives(s)
_preparation_primitives(s::Superposition) = reduce(vcat,(_preparation_primitives(t) for t in s.states))
nmodes(s::Superposition) = nmodes(first(s.states))
overlapcapability(::Superposition) = PureOverlap()
function _stateinner(s::Superposition,t::State)
    return sum(conj(c)*_stateinner(x,t) for (c,x) in zip(s.coefficients,s.states))/sqrt(BigFloat(s.normalisation))
end
function _stateinner(s::State,t::Superposition)
    return conj(_stateinner(t,s))
end
function _stateinner(s::Superposition,t::Superposition)
    return sum(conj(c)*d*_stateinner(x,y) for (c,x) in zip(s.coefficients,s.states),
               (d,y) in zip(t.coefficients,t.states))/sqrt(BigFloat(s.normalisation)*t.normalisation)
end
function logoverlap(s::Superposition,z)
    logs = [_logpoint(x,_coordinate(z)) for x in s.states]
    scale = maximum(real,logs)
    isfinite(scale) || return complex(-Inf)
    return scale+log(sum(c*exp(a-scale) for (c,a) in zip(s.coefficients,logs)))-log(s.normalisation)/2
end

# Point kernels avoid a temporary sample matrix per rejection proposal.
_qpoint!(rng,z,s::Coherent) = (z[1] = s.β+randnc(rng); z)
_qpoint!(rng,z,s::Fock) = (z[1] = sqrt(rand(rng,Gamma(s.n+1,1)))*cis(2π*rand(rng)); z)
function _qpoint!(rng,z,s::Squeezed)
    r,θ = abs(s.ϵ),angle(s.ϵ)/2
    z[1] = s.β+cis(θ)*(sqrt(exp(-r)*cosh(r)/2)*randn(rng)+im*sqrt(exp(r)*cosh(r)/2)*randn(rng))
    return z
end
function _qpoint!(rng,z,s::ProductState)
    offset = 0
    for state in s.states
        _qpoint!(rng,view(z,offset+1:offset+nmodes(state)),state)
        offset += nmodes(state)
    end
    return z
end
function _qpoint!(rng,z,s::Displaced)
    _qpoint!(rng,z,s.state)
    z .+= s.β
end
function _qpoint!(rng,z,s::Rotated)
    _qpoint!(rng,z,s.state)
    z .*= cis.(s.θ)
end
function _qpoint!(rng,z,s::Superposition)
    _qfill!(rng,reshape(z,:,1),s,SamplingDiagnostics(),_attemptlimit(1))
    return z
end
function _acceptprob(p)
    isfinite(p) && -1e-12 <= p <= 1+1e-10 || throw(ArgumentError("invalid rejection probability $p; envelope or floating-point failure."))
    return clamp(p,0,1)
end
function _qfill!(rng,Z,s::Superposition,d,limit)
    z = Vector{ComplexF64}(undef,nmodes(s))
    logs = Vector{ComplexF64}(undef,length(s.states))
    while d.accepted < size(Z,2)
        _proposal!(d,limit)
        _qpoint!(rng,z,s.states[rand(rng,s.categorical)])
        for j in eachindex(logs)
            logs[j] = _logpoint(s.states[j],z)
        end
        scale = maximum(real,logs)
        isfinite(scale) || continue
        numerator = 0.0im
        denominator = 0.0
        for j in eachindex(logs)
            A = exp(logs[j]-scale)
            numerator += s.coefficients[j]*A
            denominator += abs(s.coefficients[j])*abs2(A)
        end
        if rand(rng) < _acceptprob(abs2(numerator)/denominator)
            d.accepted += 1
            Z[:,d.accepted] .= z
        end
    end
    return Z
end

"""Finite coherent superposition of the supplied complex amplitudes (single mode)."""
CoherentSuperposition(amplitudes,coefficients) = Superposition(Coherent.(amplitudes),coefficients)
"""Four coherent components at α, iα, -α, -iα, with user-supplied coefficients."""
Compass(α;coefficients=ones(4)) = CoherentSuperposition([α,im*α,-α,-im*α],coefficients)
"""Sparse Fock sum; occupations are a vector (single mode) or modes × components integer matrix."""
FockSuperposition(n::AbstractVector{<:Integer},c) = Superposition(Fock.(n),c)
FockSuperposition(n::AbstractMatrix{<:Integer},c) = Superposition(
    Tuple(ProductState(Tuple(Fock.(col))) for col in eachcol(n)),c)
"""Two-mode `( |N,0⟩ + exp(i phase)|0,N⟩ )/√2`, N≥1; Q acceptance is 1/2."""
function NOON(N::Integer;phase=0.0)
    N>=1 || throw(ArgumentError("NOON requires N>=1."))
    isfinite(phase) || throw(ArgumentError("phase must be finite."))
    return FockSuperposition([N 0;0 N],[1,cis(phase)])
end
"""Exact `(a†)^m|β⟩` using m+1 displaced Fock components, without a coherent-state cutoff."""
function PhotonAddedCoherent(m::Integer,β::Number)
    m = _samplecount(m)
    isfinite(β) || throw(ArgumentError("β must be finite."))
    iszero(β) && return Displaced(Fock(m),0)
    logs = [loggamma(m+1)-loggamma(k+1)-loggamma(m-k+1)+(m-k)*log(abs(β))+loggamma(k+1)/2 for k in 0:m]
    c = [exp(logs[k+1]-maximum(logs))*cis(-(m-k)*angle(β)) for k in 0:m]
    return Displaced(FockSuperposition(collect(0:m),c),β)
end
"""Superposition of `D(±α)S(ϵ)|0⟩`, with the same squeezing and a relative phase."""
SqueezedCat(α,ϵ;phase=0.0) = Superposition(
    (Displaced(Squeezed(0,ϵ),α),Displaced(Squeezed(0,ϵ),-α)),[1,cis(phase)])
