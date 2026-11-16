@testset "Noise and Brownian bridge" begin
    rng = MersenneTwister(76)
    a = randnc(rng, 20_000)
    checkmean(a, 0)
    checkmean(abs2.(a), 1)
    checkmean(a.^2, 0)

    x = zeros(20_000)
    realnoise(x, nothing, 0.25, rng)
    checkmean(x, 0)
    checkmean(x.^2, 0.25)

    q, h, endpoint = 0.3, 0.5, 1.5
    realbridge(x, nothing, nothing, fill(endpoint, length(x)), q, h, rng)
    checkmean(x, q*endpoint)
    checkmean((x .- q*endpoint).^2, (1-q)*q*h)
end

@testset "Single-mode ordering and squeezing" begin
    rng = MersenneTwister(77)
    N = 45_000
    β = 0.4 + 0.7im
    r, θ, σ = 0.6, -0.92, 0.45
    for (sampler, order) in ((husimiQ, 1), (wigner, 0.5), (positiveP, 0))
        a, ap = sampler(rng, Squeezed(β, r*cis(θ)), N)
        checkmean(a, β)
        checkmean((a .- β).^2, -cis(θ)*sinh(r)*cosh(r))
        checkmean((ap .- conj(β)).*(a .- β).-order, sinh(r)^2)

        a, ap = sampler(rng, Crescent(β, r*cis(θ), σ), N)
        checkmean(a, β*exp(-σ^2/2))
        checkmean(a.^2, (β^2-cis(θ)*sinh(r)*cosh(r))*exp(-2σ^2))
        checkmean(ap.*a.-order, abs2(β)+sinh(r)^2)
    end
end

@testset "Coherent, thermal and number statistics" begin
    rng = MersenneTwister(78)
    N = 45_000
    β = 0.8 + 0.4im
    coherent = Coherent(β)
    a, ap = positiveP(coherent, 4)
    @test a == fill(β, 4)
    @test ap == conj.(a)
    @test glauberP(coherent, 4) == (a, ap)
    w, _ = wigner(rng, coherent, N)
    checkmean(w, β)
    checkmean(abs2.(w .- β), 0.5)

    # The first two radial moments distinguish circular noise from a
    # real Gaussian radius with a randomized phase.
    nbar = 0.7
    thermal = Thermal(β, nbar)
    for (sampler, offset) in ((glauberP, 0), (wigner, 0.5), (husimiQ, 1))
        a, ap = sampler(rng, thermal, N)
        @test ap == conj.(a)
        R = abs2.(a .- β)
        checkmean(R, nbar + offset)
        checkmean(R.^2, 2*(nbar + offset)^2)
    end

    @test_logs wigner(Fock(0), 10)
    @test_logs (:warn, r"Fock state sampling for W is only valid for n") wigner(Fock(1), 10)
    w, _ = wigner(rng, Fock(0), N)
    checkmean(abs2.(w), 0.5)
    checkmean(w.^2, 0)

    # Approximate Wigner branch at large occupation.
    n = 100
    w, _ = wigner(rng, Fock(n), N)
    R = abs2.(w)
    checkmean(w, 0)
    checkmean(R, n + 0.5)
    checkmean((R .- (n + 0.5)).^2, 0.25)

    # Keep the signed complex number-variance estimator.
    n = 3
    a, ap = positiveP(rng, Fock(n), N)
    R = ap.*a
    checkmean(R, n)
    checkmean(R.^2, n*(n-1))
    checkmean(R.^2 .+ (1-2n).*R .+ n^2, 0)
end

@testset "Bogoliubov ordering" begin
    rng = MersenneTwister(79)
    u, v = sqrt(2)*cis(0.4), cis(0.7)
    # Vacuum isolates the quantum offset. The thermal case checks stimulated noise.
    for (nbar, sampler, offset) in ((0, wigner, 0.5), (0, husimiQ, 1),
                                     (0, positiveP, 0), (4, glauberP, 0))
        a, ap = sampler(rng, Bogoliubov(u, v, nbar), 45_000)
        checkmean(ap.*a.-offset, 3nbar+1)
        checkmean(a.^2, (2nbar+1)*u*conj(v))
    end
    @test_throws ArgumentError glauberP(Bogoliubov(u, v, 0), 10)
end
