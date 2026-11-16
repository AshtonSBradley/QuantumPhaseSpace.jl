function fermion_annihilators(M)
    dimension = 1 << M
    return map(1:M) do j
        a = zeros(ComplexF64, dimension, dimension)
        bit = 1 << (j-1)
        for ket in 0:dimension-1
            iszero(ket & bit) && continue
            sign = iseven(count_ones(ket & (bit-1))) ? 1 : -1
            a[(ket ⊻ bit)+1, ket+1] = sign
        end
        a
    end
end

function check_fermion_moments(samples, a, density)
    (; ρ, τ) = samples
    n1 = a[1]' * a[1] - I/2
    n2 = a[2]' * a[2] - I/2
    hop = a[1]' * a[2]
    pair12 = a[1] * a[2]
    pair23 = a[2] * a[3]
    obs = (n1, n2, hermitianpart(hop), hermitianpart(-im*hop),
           hermitianpart(pair12), hermitianpart(-im*pair23))
    values = (real.(ρ[1,1,:]), real.(ρ[2,2,:]), real.(ρ[1,2,:]), imag.(ρ[1,2,:]),
              real.(τ[1,2,:]), imag.(τ[2,3,:]))
    for j in eachindex(obs)
        checkmean(values[j], real(tr(density*obs[j])))
        for i in 1:j
            sym = (obs[i]*obs[j] + obs[j]*obs[i])/2
            checkmean(values[i].*values[j], real(tr(density*sym)))
        end
    end
end

@testset "Fermionic truncated Wigner" begin
    rng = MersenneTwister(901)
    a = fermion_annihilators(3)
    vacuum = ComplexF64[1, 0, 0, 0, 0, 0, 0, 0]
    M, N = 3, 50_000

    @testset "Slater occupations and orbital basis" begin
        s = SlaterDeterminant([1, 0, 1])
        (; ρ, τ) = wigner(rng, s, N; pairing=true)
        @test all(==(0.5), ρ[1,1,:]) && all(==(-0.5), ρ[2,2,:])
        @test all(iszero, ρ[1,3,:]) && all(iszero, τ[1,2,:])
        checkmean(abs2.(ρ[1,2,:]), 0.5)
        checkmean(abs2.(τ[1,3,:]), 0.5)

        orbital = normalize(ComplexF64[1, im, 1-im])
        ket = sum(orbital[j] * (a[j]' * vacuum) for j in 1:M)
        state = SlaterDeterminant(reshape(orbital, :, 1))
        sample = wigner(rng, state, N; pairing=true)
        @test all(ishermitian, eachslice(sample.ρ; dims=3))
        @test all(t -> transpose(t) == -t, eachslice(sample.τ; dims=3))
        @test all(r -> tr(r) ≈ -0.5, eachslice(sample.ρ; dims=3))
        check_fermion_moments(sample, a, ket*ket')
    end

    @testset "Free thermal ensembles" begin
        β, μ = 1.3, 0.25
        h = [0.2 0.4im 0.1; -0.4im -0.5 0.25; 0.1 0.25 0.8]
        manybody_h = sum(h[i,j] * (a[i]' * a[j]) for i in 1:M, j in 1:M)
        number = sum(op' * op for op in a)
        gibbs = exp(-β*Hermitian(manybody_h - μ*number))
        state = FermionThermal(h, β; μ)
        check_fermion_moments(wigner(rng, state, N; pairing=true), a, gibbs/tr(gibbs))

        energies = [0.7, -0.8, 0.15]
        particles = 2
        # Brute-force enumeration includes all three configurations in the sector.
        weights = [count_ones(k) == particles ?
                   exp(-β*sum(energies[j]*Int(!iszero(k & (1 << (j-1)))) for j in 1:M)) : 0.0
                   for k in 0:7]
        state = FermionThermal(energies, β; particles)
        sample = wigner(rng, state, N; pairing=true)
        @test all(r -> tr(r) == particles-M/2, eachslice(sample.ρ; dims=3))
        check_fermion_moments(sample, a, Diagonal(weights/sum(weights)))
        @test all(==(0.5), wigner(rng, FermionThermal(energies, 1e6; particles), 100)[2,2,:])
        @test all(==(-0.5), wigner(rng, FermionThermal(energies, 1e6; particles), 100)[1,1,:])
    end

    @testset "Sampler API" begin
        s = SlaterDeterminant([0, 1])
        @test nmodes(s) == 2
        @test nmodes(FermionThermal([1, 2], 0; particles=1)) == 2
        @test size(wigner(s, 0)) == (2, 2, 0)
        @test size(wigner(rng, s, 0; pairing=true).τ) == (2, 2, 0)
        @test wigner(MersenneTwister(11), s, 10) == wigner(MersenneTwister(11), s, 10)
        @test_throws ArgumentError wigner(rng, s, -1)
        @test_throws ArgumentError SlaterDeterminant([0, 2])
        @test_throws ArgumentError SlaterDeterminant([1.0 1.0; 0.0 0.0])
        @test_throws ArgumentError FermionThermal([1.0 0.0; 1.0 0.0], 1)
        @test_throws ArgumentError FermionThermal([0, 1], 1; particles=3)
        @test_throws ArgumentError FermionThermal([0, 1], 1; particles=1, μ=0)
    end
end
