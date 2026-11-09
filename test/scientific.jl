# Compare real and imaginary sample means independently at six standard errors.
function checkmean(x, target; nsigma=6, atol=1e-10)
    @test abs(real(mean(x) - target)) <= nsigma*std(real.(x))/sqrt(length(x)) + atol
    @test abs(imag(mean(x) - target)) <= nsigma*std(imag.(x))/sqrt(length(x)) + atol
end

@testset "Scientific regressions" begin
    rng = MersenneTwister(77)
    N = 120_000
    β = 0.4 + 0.7im
    r = 0.6
    for θ in (0.43, 1.7, -0.92), (sampler, order) in ((husimiQ,1), (wigner,0.5), (positiveP,0))
        s = Squeezed(β, r*cis(θ))
        a, ap = sampler(rng,s,N)
        checkmean(a,β)
        checkmean((a.-β).^2,-cis(θ)*sinh(r)*cosh(r))
        checkmean((ap.-conj(β)).*(a.-β).-order,sinh(r)^2)
    end
    for q in (0.0,0.45), (sampler,order) in ((husimiQ,1),(wigner,0.5),(positiveP,0))
        a, ap = sampler(rng,Crescent(β,r*cis(0.8),q),N)
        checkmean(a,β*exp(-q^2/2))
        checkmean(a.^2,(β^2-cis(0.8)*sinh(r)*cosh(r))*exp(-2q^2))
        checkmean(ap.*a.-order,abs2(β)+sinh(r)^2)
    end
    for nbar in (4.0,8.0)
        u,v = sqrt(2)*cis(0.4),cis(0.7)
        a,ap = glauberP(rng,Bogoliubov(u,v,nbar),N)
        checkmean(ap.*a,3nbar+1)
        checkmean(a.^2,(2nbar+1)*u*conj(v))
    end
    @test_throws ArgumentError glauberP(Bogoliubov(sqrt(2),1,0),10)
    a,ap = wigner(rng,Fock(0),N)
    @test ap == conj.(a)
    checkmean(abs2.(a),0.5)
    checkmean(a.^2,0)
end
