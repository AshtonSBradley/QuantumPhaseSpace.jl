function _covfactor(V,tol;required=true)
    E = eigen(Symmetric(V))
    if minimum(E.values)<-tol
        required && throw(ArgumentError("covariance is not positive semidefinite."))
        return nothing
    end
    return E.vectors*Diagonal(sqrt.(max.(E.values,0)))
end
"""
    Gaussian(β, V)

Multimode Gaussian state with complex displacement vector β and real symmetric
quadrature covariance V, ordered `(x₁,p₁,x₂,p₂,...)`. Here `a=(x+ip)/√2`,
`[x,p]=i`, and `V[j,k]=⟨{δR[j],δR[k]}⟩/2`, so vacuum has `V=I/2`.
Validates `V+iΩ/2 ≥ 0` and caches W/Q/P factors. Ordinary P is available only
when `V-I/2 ≥ 0`, including semidefinite limits. Positive-P supports
`method=:drummond` for direct Gaussian sampling and `method=:canonical` for Q conversion.
"""
struct Gaussian <: State
    β::Vector{ComplexF64}
    V::Matrix{Float64}
    wfactor::Matrix{Float64}
    qfactor::Matrix{Float64}
    pfactor::Union{Nothing,Matrix{Float64}}
    function Gaussian(β::AbstractVector,V::AbstractMatrix{<:Real})
        β = ComplexF64.(collect(β))
        M = length(β)
        M>0 && size(V)==(2M,2M) || throw(DimensionMismatch("V must be 2M × 2M for M nonzero modes."))
        all(isfinite,β) && all(isfinite,V) || throw(ArgumentError("Gaussian parameters must be finite."))
        T = eltype(V)
        tol = 100*(T<:AbstractFloat ? Float64(eps(T)) : eps(Float64))*max(1,opnorm(V))
        isapprox(V,V';atol=tol,rtol=0) || throw(ArgumentError("V must be symmetric."))
        V = Matrix{Float64}((V+V')/2)
        Ω = zeros(2M,2M)
        for m in 1:M
            Ω[2m-1,2m] = 1
            Ω[2m,2m-1] = -1
        end
        eigmin(Hermitian(V+im*Ω/2))>=-tol || throw(ArgumentError("V violates V+iΩ/2 ≥ 0."))
        new(β,V,_covfactor(V,tol),_covfactor(V+I/2,tol),_covfactor(V-I/2,tol;required=false))
    end
end
Gaussian(β::Number,V::AbstractMatrix{<:Real}) = Gaussian([β],V)
nmodes(s::Gaussian) = length(s.β)
function _gaussianfill!(rng,Z,s,factor)
    M = nmodes(s)
    ξ,x = zeros(2M),zeros(2M)
    for col in eachcol(Z)
        randn!(rng,ξ)
        mul!(x,factor,ξ)
        for m in 1:M
            col[m] = s.β[m]+complex(x[2m-1],x[2m])/sqrt(2)
        end
    end
    return Z
end
function _qfill!(rng,Z,s::Gaussian,d,limit)
    _gaussianfill!(rng,Z,s,s.qfactor)
    d.proposals += size(Z,2)
    d.accepted += size(Z,2)
    return Z
end
function wigner(rng::AbstractRNG,s::Gaussian,N)
    Z = Matrix{ComplexF64}(undef,nmodes(s),_samplecount(N))
    _gaussianfill!(rng,Z,s,s.wfactor)
    return _legacy(s,Z,conj.(Z))
end
function glauberP(rng::AbstractRNG,s::Gaussian,N)
    N = _samplecount(N)
    isnothing(s.pfactor) && throw(ArgumentError("V-I/2 is not positive semidefinite; this state has no nonnegative ordinary Glauber P. Use positiveP."))
    Z = Matrix{ComplexF64}(undef,nmodes(s),N)
    _gaussianfill!(rng,Z,s,s.pfactor)
    return _legacy(s,Z,conj.(Z))
end

# The covariance of the complex quadratures is L * transpose(L).
# Negative eigenvalues contribute purely imaginary noise.
function _normal_covfactor(C)
    E = eigen(Symmetric(C))
    return E.vectors * Diagonal(sqrt.(complex.(E.values)))
end

function positiveP(rng::AbstractRNG,s::Union{Gaussian,Squeezed,SqueezedTwoMode,Bogoliubov},N;
                   method=:canonical)
    method === :canonical && return _canonical_positiveP(rng,s,N)
    method === :drummond || throw(ArgumentError("Gaussian positiveP method must be :canonical or :drummond."))

    N = _samplecount(N)
    M = nmodes(s)
    A = Matrix{ComplexF64}(undef,M,N)
    B = similar(A)
    if N > 0
        _drummond_fill!(rng,A,B,s isa Gaussian ? s : Gaussian(s))
    end
    return _legacy(s,A,B)
end

function _drummond_fill!(rng,A,B,g::Gaussian)
    M = nmodes(g)
    L = _normal_covfactor(g.V-I/2)
    ξ = Vector{Float64}(undef,2M)
    z = Vector{ComplexF64}(undef,2M)
    invsqrt2 = inv(sqrt(2.0))
    for col in axes(A,2)
        randn!(rng,ξ)
        mul!(z,L,ξ)
        for m in 1:M
            x,ip = z[2m-1],im*z[2m]
            A[m,col] = g.β[m] + invsqrt2*(x+ip)
            B[m,col] = conj(g.β[m]) + invsqrt2*(x-ip)
        end
    end
    return A,B
end

"""`D(β) S(ϵ) ρthermal(nbar) S(ϵ)† D(β)†`, using the single-mode squeezing convention."""
function SqueezedThermal(β,ϵ,nbar::Real)
    all(isfinite,(β,ϵ,nbar)) && nbar>=0 || throw(ArgumentError("finite parameters and nbar≥0 are required."))
    r,θ = abs(ϵ),angle(ϵ)/2
    R = [cos(θ) -sin(θ);sin(θ) cos(θ)]
    V = (nbar+0.5)*R*Diagonal([exp(-2r),exp(2r)])*R'
    return Gaussian([β],V)
end
"""
    GaussianBogoliubov(U,V,nbar; β=zeros(M))

Gaussian input thermal modes transformed by `a=U*b+conj(V)*b†`. `U,V` must
satisfy both bosonic canonical identities. nbar is scalar or per-input-mode.
"""
function GaussianBogoliubov(U::AbstractMatrix,V::AbstractMatrix,nbar;β=zeros(ComplexF64,size(U,1)))
    M = size(U,1)
    size(U)==size(V)==(M,M) || throw(DimensionMismatch("U,V must be square with equal dimensions."))
    U,B = ComplexF64.(U),conj.(ComplexF64.(V))
    all(isfinite,U) && all(isfinite,B) || throw(ArgumentError("U,V must be finite."))
    isapprox(U*U'-B*B',I;atol=1e-10,rtol=1e-10) &&
        isapprox(U*transpose(B),B*transpose(U);atol=1e-10,rtol=1e-10) ||
        throw(ArgumentError("U,V violate bosonic canonical identities."))
    n = _realmodevalues(nbar,M,"thermal occupation")
    S = zeros(2M,2M)
    for j in 1:M,k in 1:M
        S[2j-1,2k-1] = real(U[j,k]+B[j,k])
        S[2j-1,2k] = -imag(U[j,k]-B[j,k])
        S[2j,2k-1] = imag(U[j,k]+B[j,k])
        S[2j,2k] = real(U[j,k]-B[j,k])
    end
    return Gaussian(β,S*Diagonal(repeat(n.+0.5;inner=2))*S')
end
"""
    TwoModeSqueezedThermal(r,ϕ,n1=0,n2=n1)

Two-mode thermal squeezing with `a₁=cosh(r)b₁+exp(2iϕ)sinh(r)b₂†` and the
corresponding second mode. This preserves the legacy `SqueezedTwoMode` phase.
"""
function TwoModeSqueezedThermal(r::Real,ϕ::Real,n1::Real=0,n2::Real=n1)
    all(isfinite,(r,ϕ)) && r>=0 || throw(ArgumentError("finite r≥0 and phase are required."))
    U = Matrix{ComplexF64}(I,2,2)*cosh(r)
    B = [0 cis(2ϕ);cis(2ϕ) 0]*sinh(r)
    return GaussianBogoliubov(U,conj.(B),[n1,n2])
end
Gaussian(s::Squeezed) = SqueezedThermal(s.β,s.ϵ,0)
Gaussian(s::Thermal) = Gaussian([s.β],Matrix{Float64}(I,2,2)*(s.n̄+0.5))
Gaussian(s::Coherent) = Gaussian(Thermal(s.β,0))
Gaussian(s::Bogoliubov) = GaussianBogoliubov(reshape([s.u],1,1),reshape([s.v],1,1),s.n̄)
Gaussian(s::SqueezedTwoMode) = TwoModeSqueezedThermal(s.r,s.ϕ)
function wigner(rng::AbstractRNG,s::SqueezedTwoMode,N)
    Z = Matrix{ComplexF64}(undef,2,_samplecount(N))
    g = Gaussian(s)
    _gaussianfill!(rng,Z,g,g.wfactor)
    return _legacy(s,Z,conj.(Z))
end

"""
    PhotonAddedGaussian(state, m=1)

Normalised photon addition to a centred single-mode Gaussian. The exact Q law
is proportional to `|z|^(2m) Qstate(z)`, sampled as a finite nonnegative mixture
in the principal axes of `V+I/2`. Displaced/multimode inputs are rejected.
"""
struct PhotonAddedGaussian <: State
    state::Gaussian
    m::Int
    axes::Matrix{Float64}
    probabilities::Vector{Float64}
    categorical::Categorical{Float64,Vector{Float64}}
    radial1::Vector{Gamma{Float64}}
    radial2::Vector{Gamma{Float64}}
    function PhotonAddedGaussian(s::Gaussian,m::Integer=1)
        nmodes(s)==1 && all(iszero,s.β) || throw(ArgumentError("centred single-mode Gaussian required; displaced formula is different."))
        m = _samplecount(m)
        E = eigen(Symmetric(s.V+I/2))
        v1,v2 = E.values
        logs = [loggamma(m+1)-loggamma(k+1)-loggamma(m-k+1)+
            k*log(2v1)+(m-k)*log(2v2)+loggamma(k+0.5)+loggamma(m-k+0.5) for k in 0:m]
        p = _probabilities(exp.(logs.-maximum(logs)))
        new(s,m,E.vectors,p,Categorical(p),[Gamma(k+0.5,2v1) for k in 0:m],
            [Gamma(m-k+0.5,2v2) for k in 0:m])
    end
end
PhotonAddedGaussian(s::State,m::Integer=1) = PhotonAddedGaussian(Gaussian(s),m)
function _qfill!(rng,Z,s::PhotonAddedGaussian,d,limit)
    for col in eachcol(Z)
        j = rand(rng,s.categorical)
        x = sqrt(rand(rng,s.radial1[j]))*(rand(rng,Bool) ? 1 : -1)
        y = sqrt(rand(rng,s.radial2[j]))*(rand(rng,Bool) ? 1 : -1)
        col[1] = complex(s.axes[1,1]*x+s.axes[1,2]*y,s.axes[2,1]*x+s.axes[2,2]*y)/sqrt(2)
    end
    d.proposals += size(Z,2)
    d.accepted += size(Z,2)
    return Z
end
