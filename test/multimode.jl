@testset "Passive mixing and collective Fock states" begin
    rng = MersenneTwister(404)
    N = 50_000
    U = [1 im; im 1]/sqrt(2)
    input = ProductState(Fock(2), Fock(0))
    mixed = Passive(input, U)
    Z = sampleQ(rng, mixed, N)
    checkmean(abs2.(Z[1,:]), 2)
    checkmean(Z[1,:].*conj.(Z[2,:]), -im)
    z = [0.3im, 0.2]
    @test logoverlap(mixed, z) ≈ logoverlap(input, U'*z)

    n, u = 3, [sqrt(0.4), sqrt(0.6)*cis(0.7)]
    s = FixedNumberCondensate(n, u)
    Z = sampleQ(rng, s, N)
    @test size(Z) == (2, N)
    ℓ = abs2(u[2])
    checkmean(Z[1,:].*conj.(Z[2,:]), n*u[1]*conj(u[2]))
    checkmean(abs2.(Z[2,:]), 1+n*ℓ)
    checkmean(abs2.(Z[2,:]).^2, n*(n-1)*ℓ^2+4n*ℓ+2)

    for method in (:specialised, :generic)
        a, ap = positiveP(rng, s, N; method)
        checkmean(ap[1,:].*a[2,:], n*conj(u[1])*u[2])
        f = real.(ap[2,:].*a[2,:])
        checkmean(f, n*ℓ)
        # Check the manuscript's variance result separately from the mean.
        expected = method == :specialised ? ℓ^2*(n+2) : 2+n*ℓ*(2-ℓ)
        checkmean((f .- n*ℓ).^2, expected)
        if method == :specialised
            @test a[2,:] ≈ (u[2]/u[1]).*a[1,:]
        end
    end

    # Independent two-mode expansion from the binomial theorem.
    occupations = [n-k for k in 0:n]
    coeffs = [sqrt(binomial(n,k))*u[1]^(n-k)*u[2]^k for k in 0:n]
    expanded = FockSuperposition(vcat(transpose(occupations), transpose(collect(0:n))), coeffs)
    @test exp(logoverlap(s,z)) ≈ exp(logoverlap(expanded,z)) atol=1e-13
    vacuum = FixedNumberCondensate(0,u)
    @test all(A -> all(iszero,A), positiveP(rng,vacuum,3))
    @test size(sampleQ(rng,FixedNumberCondensate(1,[1]),3)) == (1,3)
    @test_throws ArgumentError FixedNumberCondensate(1,[1,1])
    @test_throws ArgumentError FixedNumberCondensate(-1,[1])
    @test_throws ArgumentError Superposition((mixed,),[1])
    @test_throws ArgumentError Passive(input,[1 1;0 1])
    @test_throws DimensionMismatch Passive(input,ones(1,1))
end
