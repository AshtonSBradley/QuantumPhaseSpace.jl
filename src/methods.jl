"""
    α,α⁺ = wigner(state <: State,N)

Generate `N` samples from wigner phase-space distribution for `state`.

Moments of the Wigner distribution generate symmetrically ordered quantum operator averages.
"""
function wigner(rng::AbstractRNG,state::Coherent,N)
    N = _samplecount(N)
    @unpack β = state
    α = β .+ randnc(rng,N)/sqrt(2)
    α⁺ = conj(α)
    return α, α⁺
end

"""
    α,ᾱ = husimiQ(state <: State,N)

Generate `N` samples from the Husimi-Q phase-space distribution for `state`.

Moments of the Husimi-Q distribution generate quantum operator averages that are anti-normally ordered.
"""
function husimiQ(rng::AbstractRNG,state::Coherent,N)
    N = _samplecount(N)
    @unpack β = state
    α = β .+ randnc(rng,N)
    α⁺ = conj(α)
    return α, α⁺
end

"""
    α,α⁺ = positiveP(state <: State,N)

Generate `N` samples from the positive-P (+P) phase-space distribution for `state`.

Moments of the +P distribution generate quantum operator averages that are normally ordered.

The doubled coordinates α,α⁺ need not be conjugates and are generally correlated. 
"""
function positiveP(rng::AbstractRNG,state::State,N)
    N = _samplecount(N)
    μ = sampleQ(rng,state,N)
    γ = randnc(rng,size(μ)...)
    return _legacy(state,μ .+ γ,conj.(μ .- γ))
end

function positiveP(rng::AbstractRNG,state::Coherent,N)
    N = _samplecount(N)
    @unpack β = state
    α = β*ones(N)
    α⁺ = conj(α)
    return α,α⁺
end

"""
    α,α⁺ = glauberP(state <: State,N)

Generate `N` samples from the Glauber-Sudarshan-P phase-space distribution for `state`.

Moments of the Glauber-Sudarshan-P distribution generate quantum operator averages that are normally ordered.
"""
glauberP(rng::AbstractRNG,state::Coherent,N) = positiveP(rng,state,N)

function glauberP(rng::AbstractRNG,state::Thermal,N)
    N = _samplecount(N)
    @unpack β,n̄ = state
    α = β .+ sqrt(n̄)*randnc(rng,N)
    α⁺ = conj(α)
    return α, α⁺
end

function wigner(rng::AbstractRNG,state::Thermal,N)
    N = _samplecount(N)
    @unpack β,n̄ = state
    α = β .+ sqrt(n̄ + .5)*randnc(rng,N)
    α⁺ = conj(α)
    return α, α⁺
end

function husimiQ(rng::AbstractRNG,state::Thermal,N)
    N = _samplecount(N)
    @unpack β,n̄ = state
    α = β .+ sqrt(n̄ + 1.0)*randnc(rng,N)
    α⁺ = conj(α)
    return α, α⁺
end

function husimiQ(rng::AbstractRNG,state::Squeezed,N)
    N = _samplecount(N)
    @unpack β,ϵ = state
    r = abs(ϵ)
    ϕ = angle(ϵ)/2
    ν = sqrt(exp(-r)*cosh(r)/2)*randn(rng,N) .+ im*sqrt(exp(r)*cosh(r)/2)*randn(rng,N)
    α = β .+ exp(im*ϕ)*ν
    α⁺ = conj(α)
    return α, α⁺
end

function husimiQ(rng::AbstractRNG,state::SqueezedTwoMode,N)
    N = _samplecount(N)
    @unpack r,ϕ = state
    λ = tanh(r)
    x = randnc(rng,N)*sqrt(2/(1-λ))
    y = randnc(rng,N)*sqrt(2/(1+λ))
    α = @. (-im*x +conj(y))*im*exp(im*ϕ)/2
    α⁺= conj.(α)
    β = @. (y -im*conj(x))*im*exp(im*ϕ)/2
    β⁺= conj.(β)
    return α, α⁺, β, β⁺
end

function positiveP(rng::AbstractRNG,state::SqueezedTwoMode,N)
    N = _samplecount(N)
    μa,μa⁺,μb,μb⁺ = husimiQ(rng,state,N)
    γa = randnc(rng,N)
    γb = randnc(rng,N)
    α = μa .+ γa
    α⁺ = μa⁺ .- conj.(γa)
    β = μb .+ γb
    β⁺ = μb⁺ .- conj.(γb)
    return α, α⁺, β, β⁺
end

function wigner(rng::AbstractRNG,state::Squeezed,N)
    N = _samplecount(N)
    @unpack β,ϵ = state
    r = abs(ϵ)
    ϕ = angle(ϵ)/2
    α = β .+ 0.5*(randn(rng,N)*exp(-r) .+ im*randn(rng,N)*exp(r))*exp(im*ϕ)
    α⁺ = conj(α)
    return α, α⁺
end

function husimiQ(rng::AbstractRNG,state::Crescent,N)
    N = _samplecount(N)
    α, _ = wigner(rng,state,N)
    α .+= randnc(rng,N)/sqrt(2)
    return α, conj(α)
end

function wigner(rng::AbstractRNG,state::Crescent,N)
    N = _samplecount(N)
    @unpack β,ϵ,q = state
    r = abs(ϵ)
    ϕ = angle(ϵ)/2
    α = β .+ 0.5*(randn(rng,N)*exp(-r) .+ im*randn(rng,N)*exp(r))*exp(im*ϕ)
    α = α.*exp.(im*q*randn(rng,N))
    α⁺ = conj(α)
    return α, α⁺
end

function husimiQ(rng::AbstractRNG,state::Fock,N)
    N = _samplecount(N)
    @unpack n = state
    d = Gamma(n+1,1)
    z = rand(rng,d,N)
    α = sqrt.(z).*exp.(2π*im*rand(rng,N))
    α⁺ = conj(α)
    return α, α⁺
end

function wigner(rng::AbstractRNG,state::Fock,N)
    N = _samplecount(N)
    @unpack n = state
    n == 0 && return wigner(rng,Coherent(0),N)
    n < 10 && @warn "Fock state sampling for W is only valid for n ≫ 1."
    p = 0.5*sqrt(2*n+1+2*sqrt(n^2+n))
    q = 1/(4*p)
    α = (p .+ q*randn(rng,N)).*exp.(2π*im*rand(rng,N))
    α⁺ = conj(α)
    return α, α⁺
end

"""
    α,α⁺ = positiveW(state <: State,N)

Generate `N` samples from the positive-W (+W) phase-space distribution for `state`.

Implemented states are

- `Fock(N)`

The moments of the +W distribution generate quantum operator averages that are symmmetrically ordered.
"""
function positiveW(rng::AbstractRNG,state::Fock,N)
    N = _samplecount(N)
    @unpack n = state
    if n<=320
    x1 = max(0,sqrt(n)-5); x2 = sqrt(n)+5
    (n==0||n==1) ? Pmax=0.71 : Pmax=0.6
    z = reject(rng,x->plaguerre.(x,n),[x1,x2],N,Pmax)
    μ = z.*exp.(2π*im*rand(rng,N))
    γ = randnc(rng,N)
    α = μ .+ γ
    α⁺ = conj(μ .- γ)
    return α, α⁺
    else
    x1 = max(0,sqrt(n)-5); x2 = sqrt(n)+5
    z = reject(rng,x->plaguerre_asymptotic.(x,n),[x1,x2],N,0.6)
    μ = z.*exp.(2π*im*rand(rng,N))
    γ = randnc(rng,N)
    α = μ .+ γ
    α⁺ = conj(μ .- γ)
    return α, α⁺
end
end

function glauberP(rng::AbstractRNG,state::Bogoliubov,N)
    return glauberP(rng,Gaussian(state),N)
end

function wigner(rng::AbstractRNG,state::Bogoliubov,N)
    N = _samplecount(N)
    @unpack u,v,n̄ = state
    β,β⁺ = wigner(rng,Thermal(0.0,n̄),N)
    α = u*β + conj(v)*β⁺
    α⁺ = conj.(α)
    return α,α⁺
end

function husimiQ(rng::AbstractRNG,state::Bogoliubov,N)
    N = _samplecount(N)
    @unpack u,v,n̄ = state
    β,β⁺ = wigner(rng,Thermal(0.0,n̄),N)
    α = u.*β .+ conj(v).*β⁺ .+ randnc(rng,N)/sqrt(2)  # convolve with Gaussian 
    α⁺ = conj.(α)
    return α,α⁺
end

# Backward-compatible default-RNG entry points.
for sampler in (:husimiQ,:wigner,:glauberP,:positiveP,:positiveW)
    @eval $sampler(state::State,N; kwargs...) = $sampler(Random.default_rng(),state,N; kwargs...)
end

positiveP(rng::AbstractRNG,state::Thermal,N) = glauberP(rng,state,N)
