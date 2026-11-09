function _coordinatematrices(s,result,N)
    if length(result)==4
        return vcat(transpose(result[1]),transpose(result[3])),
            vcat(transpose(result[2]),transpose(result[4]))
    end
    A,B = result
    return reshape(A,nmodes(s),N),reshape(B,nmodes(s),N)
end
"""
    positiveW([rng], state, N)

Exact doubled positive-W sampling. Generic states add independent circular noise
η=ξ/√2 to a valid positive-P pair as `(α+η,αplus+conj(η))`. Coordinates are
generally correlated and need not be conjugates. They are not ordinary W samples
and must not be projected into ordinary truncated-Wigner initial data.

Fock states use an exact binomial–Gamma radial mixture with no spatial window
or asymptotic switch. `method=:convolution` selects the generic construction
for a benchmark of the same representation and physical estimators.
"""
function positiveW(rng::AbstractRNG,s::State,N;method=:convolution)
    method==:convolution || throw(ArgumentError("generic positiveW supports method=:convolution."))
    N = _samplecount(N)
    A,B = _coordinatematrices(s,positiveP(rng,s,N),N)
    η = randnc(rng,nmodes(s),N)/sqrt(2)
    return _legacy(s,A.+η,B.+conj.(η))
end
function positiveW(rng::AbstractRNG,s::Fock,N;method=:direct)
    method==:convolution && return invoke(positiveW,Tuple{AbstractRNG,State,Any},rng,s,N;method)
    method==:direct || throw(ArgumentError("Fock positiveW method must be :direct or :convolution."))
    N = _samplecount(N)
    a,ap = Vector{ComplexF64}(undef,N),Vector{ComplexF64}(undef,N)
    binomial = Binomial(s.n,2/3)
    for i in eachindex(a)
        k = rand(rng,binomial)
        μ = sqrt(rand(rng,Gamma(k+1,1.5)))*cis(2π*rand(rng))
        γ = randnc(rng)
        a[i],ap[i] = μ+γ,conj(μ-γ)
    end
    return a,ap
end
