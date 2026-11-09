#helper functions
"""
    a = randnc(N...)

Return an array of dimension `length(N)`, containing samples of complex random
variables with mean zero and variance one:
```math
⟨|a|^2⟩ = 1; ⟨a^2⟩ = ⟨(a^*) ^2⟩ = 0,
```
Useful for creating a range of noises that show up in phase-space simulations.

Usage:

`a = randnc(10)` returns a `10-element Array{Complex{Float64},1}`.

`a = randnc(50,100)` returns a `50x100-element Array{Complex{Float64},2}`.
"""
randnc(rng::AbstractRNG,args...) = randn(rng,ComplexF64,args...)
randnc(args...) = randnc(Random.default_rng(),args...)


"""
    x = reject(P,w,N,Pmax)

Generate `x` distributed according to a probability distribution by rejection sampling over a finite window.

`P`: probability distribution, normalized to 1.

`w=[w1,w2]`: window for sampling `x`.

`N`: number of samples.

`Pmax`: numerical upper bound for `P`: `Pmax ≧ max(P(w))`.

"""
function reject(rng::AbstractRNG,P,w,N,Pmax;
                max_attempts=_attemptlimit(N), diagnostics=false)
    N = _samplecount(N)
    length(w) == 2 || throw(ArgumentError("reject requires a two-element window."))
    a,b = w
    isfinite(a) && isfinite(b) && a < b || throw(ArgumentError("reject requires finite w[1] < w[2]."))
    isfinite(Pmax) && Pmax > 0 || throw(ArgumentError("reject requires finite Pmax > 0."))
    max_attempts = _samplecount(max_attempts)
    samples = Vector{Float64}(undef,N)
    stats = SamplingDiagnostics()
    while stats.accepted < N
        _proposal!(stats,max_attempts)
        y = a + rand(rng)*(b-a)
        py = P(y)
        isfinite(py) && 0 <= py <= Pmax || throw(ArgumentError("reject requires finite 0 <= P(y) <= Pmax."))
        if rand(rng)*Pmax < py
            stats.accepted += 1
            samples[stats.accepted] = y
        end
    end
    return diagnostics ? (samples,stats) : samples
end
reject(P,w,N,Pmax;kwargs...) = reject(Random.default_rng(),P,w,N,Pmax;kwargs...)

# Laguerre methods for +W
plaguerre(x,n) = (4/3)*exp(-2*x^2/3 + log(x) + log(laguerren(-4*x^2/3,n)) - n*log(3))
plaguerre_asymptotic(x,n) = exp(-(x-sqrt(n+1))^2)/sqrt(pi)
laguerren(x,n) = sf_laguerre_n(n,0.0,x)


"""
Complex variable SDE's in `DifferentialEquations.jl`
have complex noises by default dispatch on complex fields.
This method defines a real noise for use with the library.
"""
function realnoise(rand_vec,W,dt,rng)
for i in eachindex(rand_vec)
    rand_vec[i] = randn(rng)
    end
    rand_vec .*= sqrt(abs(dt))
end

"""
Complex variable SDE's in `DifferentialEquations.jl`
have complex noises by default dispatch on complex fields.
This method defines a real brownian bridge for use with the library.
"""
function realbridge(rand_vec,W,W0,Wh,q,h,rng)
for i in eachindex(rand_vec)
    rand_vec[i] = randn(rng)
    end
    rand_vec .= sqrt((1 .- q).*q.*abs(h)).*rand_vec.+q.*Wh
end
