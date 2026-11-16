@testset "Finite superpositions and cats" begin
    rng = MersenneTwister(301)
    K = 48
    a = annihilator(K)
    α,φ = 0.8+0.3im,0.7
    v = coherentket(α,K)+cis(φ)*coherentket(-α,K)
    normalize!(v)
    s = CoherentSuperposition([α,-α],[1,cis(φ)])
    checkqket(rng,s,v)
    t = CoherentSuperposition([-α,α],[2im*cis(φ),2im])
    @test t.expected_acceptance ≈ s.expected_acceptance
    @test abs2(exp(logoverlap(s,0.8-0.4im))) ≈ abs2(exp(logoverlap(t,0.8-0.4im)))
    @test_throws ArgumentError Superposition((Thermal(0,1),),[1])
    @test_throws ArgumentError Superposition((Fock(0),Fock(0)),[1,-1])
    @test_throws DimensionMismatch Superposition((Fock(0),ProductState(Fock(0),Fock(0))),[1,1])
    @test_throws SamplingError sampleQ(s,10;max_attempts=0)
    for t in (s, NOON(20), Compass(α))
        Z,d = sampleQ(rng,t,15_000;diagnostics=true)
        p = t.expected_acceptance
        @test abs(acceptance_rate(d)-p)<6sqrt(p*(1-p)/d.proposals)
        @test size(Z) == (nmodes(t),15_000)
    end
    # Sparse bosonic-code-like finite sum, with interference between |0> and |2>.
    v = zeros(ComplexF64,K); v[[1,3,6]] .= [1,im,0.3]; normalize!(v)
    checkqket(rng,FockSuperposition([0,2,5],[1,im,0.3]),v)
    β = 0.4+0.5im
    v = adjoint(a)^3*coherentket(β,K)
    normalize!(v)
    checkqket(rng,PhotonAddedCoherent(3,β),v)
    @test_throws ArgumentError OddCat(0)

    # Exercise the one-photon limit and a finite odd cat where Gamma
    # rejection has a measurable effect. The large odd cat uses Gaussians.
    Aodd = 0.6^2
    for (amplitude, phase, expected_rate) in ((1e-6, 0, 1), (1e-6+2e-6im, π, 1),
                                             (0.6cis(0.4), π, (1-Aodd/3)^2*sinh(Aodd)/Aodd),
                                             (1.8im, π, -expm1(-2*1.8^2)/2),
                                             (α, φ, (1+exp(-2abs2(α))*cos(φ))/2))
        t = Cat(amplitude,phase)
        # All chosen cases accept about half the proposals or more.
        Z,d = sampleQ(rng,t,35_000;diagnostics=true,max_attempts=140_000)
        z = vec(Z)
        A,δ = abs2(amplitude),exp(-2abs2(amplitude))
        occupation = phase == π ? A/tanh(A) : A*(1-δ*cos(phase))/(1+δ*cos(phase))
        first_moment = phase == φ ? -im*amplitude*δ*sin(phase)/(1+δ*cos(phase)) : 0im
        checkmean(z,first_moment)
        checkmean(z.^2,amplitude^2)
        checkmean(abs2.(z),1+occupation)
        @test abs(acceptance_rate(d)-expected_rate) <=
              6sqrt(expected_rate*(1-expected_rate)/d.proposals)+2/d.proposals
    end
    # Complex displaced squeezed overlaps: compare with matrix exponentials.
    ϵ = 0.35cis(0.8)
    S = exp((conj(ϵ)*a^2-ϵ*adjoint(a)^2)/2)
    D = exp(β*adjoint(a)-conj(β)*a)
    v = D*S[:,1]
    z = 0.7-0.4im
    @test exp(logoverlap(Squeezed(β,ϵ),z)) ≈ dot(coherentket(z,K),v) atol=1e-12
    t = SqueezedCat(β,ϵ;phase=0.4)
    vm = exp(-β*adjoint(a)+conj(β)*a)*S[:,1]
    vs = v+cis(0.4)*vm; normalize!(vs)
    checkqket(rng,t,vs)
    @test t.normalisation ≈ norm(v+cis(0.4)*vm)^2/4 atol=1e-12
    # Displaced Fock Gram and phase, independently from exp(D).
    states = (Displaced(Fock(2),β),Displaced(Fock(3),-0.3im),Squeezed(β,ϵ))
    v1 = D[:,3]; v2 = exp(-0.3im*adjoint(a)-0.3im*a)[:,4]
    t = Superposition(states,[1,0.3im,-0.4])
    cols = hcat(v1,v2,v)
    @test t.gram ≈ cols'*cols atol=2e-12
    vs = cols*[1,0.3im,-0.4]; normalize!(vs)
    checkqket(rng,t,vs)
    Z = sampleQ(rng,NOON(2;phase=φ),55_000)
    checkmean(conj.(Z[1,:]).^2 .* Z[2,:].^2,cis(φ))
    @test real(logoverlap(Fock(10000),sqrt(10000))) > -10
    @test isfinite(real(logoverlap(CoherentSuperposition([40,-40],[1,1]),40)))
end
