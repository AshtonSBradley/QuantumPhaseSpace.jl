using QuantumPhaseSpace, Random, Statistics, LinearAlgebra, Printf

# Base.@timed, warmed callbacks, median of five runs. Sampling time includes
# output allocation; variance and estimator evaluation are outside the timer.
function generic_positive(rng,s,N)
    μ,d = sampleQ(rng,s,N;diagnostics=true)
    γ = randn(rng,ComplexF64,size(μ))
    return μ.+γ,conj.(μ.-γ),acceptance_rate(d)
end
function samplecase(rng,s,N,representation,method)
    if representation==:Q
        Z,d = sampleQ(rng,s,N;diagnostics=true)
        return Z,conj.(Z),acceptance_rate(d)
    elseif representation==:P
        if method==:generic
            A,B,acceptance = generic_positive(rng,s,N)
        else
            A,B = s isa FixedNumberCondensate ? positiveP(rng,s,N;method=:specialised) : positiveP(rng,s,N)
            acceptance = 1.0
        end
        return reshape(A,nmodes(s),N),reshape(B,nmodes(s),N),acceptance
    end
    A,B = positiveW(rng,s,N;method)
    return reshape(A,nmodes(s),N),reshape(B,nmodes(s),N),1.0
end
function benchcase(io,name,construct,representation,method,target)
    construct() # warm preparation before measuring
    preparation = [@timed construct() for _ in 1:5]
    prep_seconds = median(t.time for t in preparation)
    prep_bytes = median(t.bytes for t in preparation)
    s = construct()
    samplecase(MersenneTwister(1),s,128,representation,method)
    for N in (100,10_000,100_000)
        timings = [@timed samplecase(MersenneTwister(101+i),s,N,representation,method) for i in 1:5]
        elapsed = median(t.time for t in timings)
        bytes = median(t.bytes for t in timings)
        A,B,acceptance = last(timings).value
        order = representation==:Q ? 1.0 : representation==:Wplus ? 0.5 : 0.0
        f = A[1,:].*B[1,:].-order
        vr,vi = var(real.(f)),var(imag.(f))
        μ = mean(f)
        @printf(io,"%s\t%s\t%s\t%d\t%.9g\t%.0f\t%.9g\t%.0f\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\n",
            name,representation,method,N,prep_seconds,prep_bytes,elapsed,bytes,
            acceptance,real(μ),imag(μ),target,vr,vi,sqrt(vr/N),vr*elapsed/N)
    end
end
function run_sampling_benchmarks(path=joinpath(@__DIR__,"results.tsv"))
    u = fill(1/sqrt(16),16).*cis.(range(0,2π;length=16))
    cases = [
        ("even_cat",()->EvenCat(1.2+0.4im),:Q,:default,1.6*tanh(1.6)),
        ("even_cat",()->EvenCat(1.2+0.4im),:P,:generic,1.6*tanh(1.6)),
        ("small_odd_cat",()->OddCat(0.0001im),:Q,:default,1e-8/tanh(1e-8)),
        ("small_odd_cat",()->OddCat(0.0001im),:P,:generic,1e-8/tanh(1e-8)),
        ("added_thermal",()->PhotonAddedThermal(2,0.7),:Q,:default,4.1),
        ("added_thermal",()->PhotonAddedThermal(2,0.7),:P,:generic,4.1),
        ("added_coherent",()->PhotonAddedCoherent(1,0.8+0.3im),:Q,:default,(0.73^2+3*0.73+1)/1.73),
        ("added_coherent",()->PhotonAddedCoherent(1,0.8+0.3im),:P,:generic,(0.73^2+3*0.73+1)/1.73),
        ("NOON_2",()->NOON(2;phase=0.7),:Q,:default,1.0),
        ("NOON_2",()->NOON(2;phase=0.7),:P,:generic,1.0),
        ("condensate_30_M16",()->FixedNumberCondensate(30,u),:P,:specialised,30/16),
        ("condensate_30_M16",()->FixedNumberCondensate(30,u),:P,:generic,30/16),
        ("thermal",()->Thermal(0.4im,3),:P,:specialised,3.16),
        ("thermal",()->Thermal(0.4im,3),:P,:generic,3.16)]
    open(path,"w") do io
        println(io,"state\trepresentation\tmethod\tN\tprep_seconds\tprep_bytes\tseconds\tbytes\tacceptance\tmean_real\tmean_imag\ttarget\tvariance_real\tvariance_imag\tstderr_real\tvariance_time_per_sample")
        for case in cases
            benchcase(io,case...)
        end
        for n in (0,1,320,10000),method in (:direct,:convolution)
            benchcase(io,"Fock_$n",()->Fock(n),:Wplus,method,Float64(n))
        end
    end
    println("Benchmark results: ",path)
    return path
end
if abspath(PROGRAM_FILE)==@__FILE__
    run_sampling_benchmarks()
end
