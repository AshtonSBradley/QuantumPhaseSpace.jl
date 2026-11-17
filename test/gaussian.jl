@testset "Gaussian covariances and quadrature conventions" begin
    rng = MersenneTwister(606)
    # [x,p]=i and [X,Y]=2i describe the same a, with X=√2x,Y=√2p.
    aop = annihilator(8)
    xop,pop = (aop+aop')/sqrt(2),(aop-aop')/(im*sqrt(2))
    Xop,Yop = sqrt(2)*xop,sqrt(2)*pop
    @test (xop+im*pop)/sqrt(2) ≈ aop
    @test (Xop+im*Yop)/2 ≈ aop
    @test (xop*pop-pop*xop)[1:7,1:7] ≈ im*Matrix(I,7,7)
    @test (Xop*Yop-Yop*Xop)[1:7,1:7] ≈ 2im*Matrix(I,7,7)
    z,_ = wigner(rng,Gaussian(Coherent(0)),40_000)
    checkmean((sqrt(2)*real.(z)).^2,0.5)
    @test_throws ArgumentError Gaussian(0,[0.1 0;0 0.1])
    @test_throws ArgumentError Gaussian(0,[1.0 0.2;0 1.0])
    @test_throws ArgumentError Gaussian(0,[Inf 0;0 1])
    @test_throws DimensionMismatch Gaussian([0,0],Matrix(I,2,2))
    classical = Gaussian(0,[0.5 0;0 1.2])
    z,_ = glauberP(rng,classical,30_000)
    @test all(iszero,real.(z))
    checkmean(abs2.(z),0.35)
    @test_throws ArgumentError glauberP(Gaussian(Squeezed(0,0.3im)),1)
    r, φ = 0.4, 0.7
    N = 45_000
    s = TwoModeSqueezedThermal(r,φ,0.3,0.6)
    occupation1 = 0.3 + 1.9*sinh(r)^2
    occupation2 = 0.6 + 1.9*sinh(r)^2
    pair = 1.9*cis(2φ)*sinh(r)*cosh(r)
    for (sampler, offset) in ((wigner,0.5),(husimiQ,1),(positiveP,0))
        a, ap = sampler(rng,s,N)
        checkmean(ap[1,:].*a[1,:].-offset,occupation1)
        checkmean(ap[2,:].*a[2,:].-offset,occupation2)
        checkmean(a[1,:].*a[2,:],pair)
    end

    legacy = SqueezedTwoMode(r,φ)
    for (sampler, offset) in ((wigner,0.5),(husimiQ,1),(positiveP,0))
        a, ap, b, bp = sampler(rng,legacy,N)
        checkmean(ap.*a.-offset,sinh(r)^2)
        checkmean(a.*b,cis(2φ)*sinh(r)*cosh(r))
    end

    # An oriented single-mode covariance tests displacement and quadrature units.
    β, θ, nbar = 0.3im, 0.8, 0.7
    s = SqueezedThermal(β,r*cis(θ),nbar)
    z,_ = husimiQ(rng,s,N)
    checkmean(z,β)
    checkmean((z.-β).^2,-(2nbar+1)*cis(θ)*sinh(r)*cosh(r))
    checkmean(abs2.(z.-β),nbar+(2nbar+1)*sinh(r)^2+1)
    @test Gaussian(s.β,(2s.V)/2).V == s.V

    u,v = sqrt(2)*cis(0.2),cis(0.7)
    s = Gaussian(Bogoliubov(u,v,4))
    z,_ = wigner(rng,s,N)
    checkmean(abs2.(z),13.5)
    checkmean(z.^2,9u*conj(v))
    @test_throws ArgumentError GaussianBogoliubov(ones(2,2),zeros(2,2),0)
end

function added_gaussian_density(K,ϵ,nbar,m)
    a = annihilator(K)
    S = exp((conj(ϵ)*a^2-ϵ*adjoint(a)^2)/2)
    p = [nbar^n/(nbar+1)^(n+1) for n in 0:K-1]
    ρ = S*Diagonal(p)*S'
    ρ = adjoint(a)^m*ρ*a^m
    return ρ/tr(ρ),a
end
@testset "Photon-added centred Gaussian" begin
    rng = MersenneTwister(607)
    ϵ = 0.35cis(0.7)
    for (nbar,m) in ((0.0,0),(0.0,1),(0.4,2))
        g = SqueezedThermal(0,ϵ,nbar)
        s = PhotonAddedGaussian(g,m)
        ρ,a = added_gaussian_density(96,ϵ,nbar,m)
        ρsmall,asmall = added_gaussian_density(64,ϵ,nbar,m)
        @test abs(tr(ρ*a^2*adjoint(a)^2)-tr(ρsmall*asmall^2*adjoint(asmall)^2))<1e-7
        z = vec(sampleQ(rng,s,45_000))
        for (p,q) in ((0,2),(1,1),(2,2))
            checkmean(conj.(z).^p .* z.^q,tr(ρ*a^q*adjoint(a)^p))
        end
        if m==1
            v = eigvals(Symmetric(g.V+I/2))
            @test s.probabilities ≈ reverse(v)/sum(v)
        end
    end
    @test_throws ArgumentError PhotonAddedGaussian(Gaussian(Coherent(0.1)),1)
    @test_throws ArgumentError PhotonAddedGaussian(TwoModeSqueezedThermal(0.2,0),1)
    @test_throws ArgumentError PhotonAddedGaussian(Gaussian(Coherent(0)),-1)
end

@testset "Direct Gaussian positive-P" begin
    rng = MersenneTwister(608)
    N = 20_000
    vacuum = Gaussian(Coherent(0))
    thermal = Gaussian(Thermal(0.4+0.2im,1.3))
    legacy = SqueezedTwoMode(0.5,π/4)

    @testset "Normal covariance factor" begin
        semidefinite = [0.5 0.5; 0.5 0.5]
        nearzero = Diagonal([-eps(Float64),eps(Float64)])
        strongly_squeezed = SqueezedThermal(0,3cis(0.7),0).V-I/2
        for C in (vacuum.V-I/2,semidefinite,nearzero,strongly_squeezed)
            L = QuantumPhaseSpace._normal_covfactor(C)
            @test L*transpose(L) ≈ C
        end
    end

    @testset "Layouts, RNG and classical limits" begin
        @test all(isempty,positiveP(rng,vacuum,0;method=:drummond))
        @test map(length,positiveP(rng,legacy,0;method=:drummond)) == (0,0,0,0)
        @test map(size,positiveP(rng,Gaussian(legacy),0;method=:drummond)) == ((2,0),(2,0))
        @test_throws ArgumentError positiveP(rng,vacuum,-1;method=:drummond)
        @test_throws ArgumentError positiveP(rng,legacy,1;method=:unknown)
        @test positiveP(MersenneTwister(1),legacy,16;method=:drummond) ==
              positiveP(MersenneTwister(1),legacy,16;method=:drummond)
        for s in (legacy,Gaussian(legacy))
            @test positiveP(MersenneTwister(2),s,16;method=:canonical) ==
                  positiveP(MersenneTwister(2),s,16)
        end
        bog = Bogoliubov(sqrt(2),cis(0.5),0.2)
        @test positiveP(MersenneTwister(3),bog,16;method=:drummond) ==
              positiveP(MersenneTwister(3),Gaussian(bog),16;method=:drummond)
        for β in (0,0.2+0.3im)
            @test positiveP(Gaussian(Coherent(β)),16;method=:drummond) ==
                  (fill(β,16),fill(conj(β),16))
        end
        α,αp = positiveP(rng,thermal,N;method=:drummond)
        @test αp == conj.(α)
        n = abs2.(α.-only(thermal.β))
        checkmean(n,1.3)
        checkmean(n.^2,2*1.3^2)
    end

    @testset "Displaced single-mode squeezing" begin
        for (r,θ,β) in ((0.15,0.0,0.0),(0.8,π/3,0.6+0.2im),(2.0,-π/2,-0.4im))
            a,ap = positiveP(rng,Squeezed(β,r*cis(θ)),N;method=:drummond)
            n = sinh(r)^2
            m = -cis(θ)*sinh(r)*cosh(r)
            checkmean(a,β)
            checkmean(ap.*a,n+abs2(β))
            checkmean(a.^2,m+β^2)
            checkmean((ap.*a).^2,2n^2+abs2(m)+4n*abs2(β)+2real(conj(β)^2*m)+abs2(β)^2)
        end
    end

    @testset "Two-mode squeezing" begin
        for (r,ϕ) in ((0.4,0.0),(0.8,π/2),(2.0,π/4))
            a,ap,b,bp = positiveP(rng,SqueezedTwoMode(r,ϕ),N;method=:drummond)
            n = sinh(r)^2
            m = cis(2ϕ)*sinh(r)*cosh(r)
            na,nb = ap.*a,bp.*b
            checkmean(na,n)
            checkmean(nb,n)
            checkmean(a.*b,m)
            # The strongly squeezed complex-phase case also tests Wick
            # contractions and minimum nonconjugacy.
            if r == 2
                checkmean(na.^2,2n^2)
                checkmean(na.*nb,n^2+abs2(m))
                checkmean((a.*b).^2,2m^2)
                checkmean(abs2.(a.-conj.(ap)),-expm1(-2r))
                checkmean(abs2.(b.-conj.(bp)),-expm1(-2r))
            end
        end

        r,d = 0.65,1.4
        g = Gaussian(fill(d,2),TwoModeSqueezedThermal(r,0).V)
        a,ap = positiveP(rng,g,N;method=:drummond)
        na,nb = ap[1,:].*a[1,:],ap[2,:].*a[2,:]
        Δ = na.-nb
        # The mean difference is zero. Restore the commutator in (n₁−n₂)².
        checkmean(Δ.^2 .+ na .+ nb,2d^2*exp(-2r))
    end

    @testset "Correlated multimode squeezing" begin
        M = 3
        mix = [cis(2π*j*k/M)/sqrt(M) for j in 0:M-1,k in 0:M-1]
        r = [0.2,0.65,0.85]
        θ = [0.1,-0.8,1.1]
        nthermal = [0.0,0.3,0.15]
        β = [0.25im,-0.4,0.3-0.2im]
        U = mix*Diagonal(cosh.(r))
        V = conj.(mix*Diagonal(-cis.(θ).*sinh.(r)))
        g = GaussianBogoliubov(U,V,nthermal;β)
        A,B = positiveP(rng,g,N;method=:drummond)
        @test size(A) == size(B) == (M,N)

        n = nthermal.+(2nthermal.+1).*sinh.(r).^2
        m = -(2nthermal.+1).*cis.(θ).*sinh.(r).*cosh.(r)
        normal = conj.(mix)*Diagonal(n)*transpose(mix)
        anomalous = mix*Diagonal(m)*transpose(mix)
        for i in 1:M
            checkmean(A[i,:],β[i])
            checkmean(B[i,:],conj(β[i]))
        end
        A .-= β
        B .-= conj.(β)
        for i in 1:M,j in 1:M
            checkmean(B[i,:].*A[j,:],normal[i,j])
        end
        # The anomalous products are symmetric trajectory by trajectory.
        for i in 1:M,j in i:M
            checkmean(A[i,:].*A[j,:],anomalous[i,j])
            checkmean(B[i,:].*B[j,:],conj(anomalous[i,j]))
        end
        checkmean(B[1,:].*B[2,:].*A[1,:].*A[2,:],
                  normal[1,1]*normal[2,2]+abs2(normal[1,2])+abs2(anomalous[1,2]))
        checkmean(A[1,:].^3 .* A[3,:],3anomalous[1,1]*anomalous[1,3])
    end
end
