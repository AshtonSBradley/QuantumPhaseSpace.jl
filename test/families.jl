@testset "Exact Q families" begin
    rng = MersenneTwister(203)
    N = 100_000
    β = 0.7+0.2im
    a,_ = husimiQ(rng,Displaced(Fock(3),β),N)
    checkmean(a,β)
    checkmean(abs2.(a.-β),4)
    s = Displaced(Coherent(0.2im),β)
    z = 0.4-0.7im
    phase = (β*conj(z)-conj(β)*z)/2
    @test exp(logoverlap(s,z)) ≈ exp(phase+logoverlap(s.state,z-β))
    for m in (0,1,3), nbar in (0.0,0.7)
        s = PhotonAddedThermal(m,nbar)
        a,_ = husimiQ(rng,s,N)
        checkmean(abs2.(a),(m+1)*(nbar+1))
        checkmean(abs2.(a).^2,(m+1)*(m+2)*(nbar+1)^2)
        a,ap = positiveP(rng,s,N)
        checkmean(ap.*a,(m+1)*nbar+m)
    end
    for m in (0,1,3)
        s = PhotonSubtractedThermal(m,0.7)
        for (sampler,order) in ((glauberP,0),(positiveP,0),(wigner,0.5),(husimiQ,1))
            a,ap = sampler(rng,s,N)
            @test ap == conj.(a)
            checkmean(ap.*a.-order,(m+1)*0.7)
        end
        a,ap = positiveP(rng,s,N)
        checkmean(ap.^2 .* a.^2,(m+1)*(m+2)*0.7^2)
    end
    @test_throws ArgumentError PhotonSubtractedThermal(1,0)
    @test_throws ArgumentError PhotonAddedThermal(-1,1)
    @test_throws ArgumentError PhotonAddedThermal(0,Inf)
    @test_throws ArgumentError Mixture([Fock(0)],[0])
    @test_throws DimensionMismatch Mixture([Fock(0),ProductState(Fock(1),Fock(0))],[1,1])
    p = ProductState(Coherent(β),Fock(2))
    Z = sampleQ(rng,p,N)
    @test size(Z) == (2,N)
    checkmean(Z[1,:],β)
    checkmean(abs2.(Z[2,:]),3)
    @test logoverlap(p,[z,z]) ≈ logoverlap(Coherent(β),z)+logoverlap(Fock(2),z)
    a,_ = husimiQ(rng,NumberMixture([0.25,0,0.75]),N)
    checkmean(abs2.(a),2.5)
    a,_ = husimiQ(rng,Rotated(Coherent(β),0.7),N)
    checkmean(a,β*cis(0.7))
    for s in (p,Displaced(Fock(0),β),Rotated(Fock(1),0.2),NumberMixture([1]),
              PhotonAddedThermal(0,0),PhotonSubtractedThermal(0,0))
        @test size(sampleQ(rng,s,0)) == (nmodes(s),0)
    end
end
