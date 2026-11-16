using QuantumPhaseSpace, Random, Statistics

"""Run the exact-state examples with a reproducible RNG and return the samples."""
function exact_sampling_examples(;N=10_000,rng=MersenneTwister(20261110))
    cat = EvenCat(1.2+0.4im)
    added_thermal = PhotonAddedThermal(2,0.7)
    added_coherent = PhotonAddedCoherent(1,0.8+0.3im)
    noon = NOON(2;phase=π/3)
    u = [sqrt(0.3),sqrt(0.7)*cis(0.4)]
    condensate = FixedNumberCondensate(20,u)
    lossy_cat = Loss(cat,0.65)
    diffused_squeezed = PhaseDiffused(Squeezed(0.5im,0.4cis(0.8)),0.3)
    added_squeezed = PhotonAddedGaussian(Squeezed(0,0.4cis(0.7)))
    displaced_fock = Displaced(Fock(3),0.5im)
    rotated = Rotated(Coherent(0.4+0.2im),0.7)
    number_mixture = NumberMixture([0.1,0.3,0.6])
    classical_mixture = Mixture((Thermal(0,0.5),Coherent(1im)),[0.6,0.4])
    product = ProductState(Fock(2),Coherent(0.5))
    subtracted_thermal = PhotonSubtractedThermal(1,0.7)
    compass = Compass(0.7+0.2im)
    sparse_code = FockSuperposition([0,2,4],[1,0.5im,-0.3])
    coherent_sum = CoherentSuperposition([0.5,0.5im,-0.5],[1,im,0.4])
    squeezed_cat = SqueezedCat(0.5+0.2im,0.3cis(0.7);phase=0.4)
    passive = Passive(ProductState(Fock(2),Fock(0)),[1 im;im 1]/sqrt(2))
    squeezed_thermal = SqueezedThermal(0.2im,0.3cis(0.7),0.7)
    two_mode_thermal = TwoModeSqueezedThermal(0.3,0.7,0.5,0.7)
    bogoliubov = GaussianBogoliubov(cosh(0.3)*[1 0;0 1],
        sinh(0.3)*[0 cis(-0.7);cis(-0.7) 0],[0.2,0.4])
    covariance_state = Gaussian([0.2+0.1im],[0.8 0.1;0.1 0.9])
    # Every Q result is modes × samples, including a single-mode row.
    states = (;cat,added_thermal,added_coherent,noon,condensate,lossy_cat,
              diffused_squeezed,added_squeezed,displaced_fock,rotated,
              number_mixture,classical_mixture,product,subtracted_thermal,
              compass,sparse_code,coherent_sum,squeezed_cat,passive,
              squeezed_thermal,two_mode_thermal,bogoliubov,covariance_state)
    Q = map(s->sampleQ(rng,s,N),states)
    # Single-mode positive-P returns vectors; new multimode states return matrices.
    P = map(s->positiveP(rng,s,N),states)
    _,diagnostics = sampleQ(rng,cat,N;diagnostics=true)
    println("Even-cat Q acceptance: ",acceptance_rate(diagnostics))
    println("Lossy-cat occupation from Q: ",mean(abs2,Q.lossy_cat)-1)
    return (;states,Q,P)
end

if abspath(PROGRAM_FILE)==@__FILE__
    exact_sampling_examples()
end
