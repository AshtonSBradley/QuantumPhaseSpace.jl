using QuantumPhaseSpace, Random, Statistics, Printf

standard_errors(z) = (std(real.(z)),std(imag.(z))) ./ sqrt(length(z))

function benchmark_two_mode(N=200_000)
    println("N = $N, Julia $VERSION, ", Sys.cpu_info()[1].model, ", Julia threads = ", Threads.nthreads())
    println("r\tmethod\tD1\tD2\tSE Re n1\tSE Im n1\tSE Re n2\tSE Im n2\tSE Re ab\tSE Im ab\tSE Re Var diff\tSE Im Var diff\tms\tbytes/sample")
    for r in (0.0,0.5,1.0,2.0,3.0), method in (:canonical,:drummond)
        s = SqueezedTwoMode(r,0)
        positiveP(Xoshiro(1),s,N;method)
        timings = map(1:9) do seed
            rng = Xoshiro(seed)
            GC.gc()
            run = @timed positiveP(rng,s,N;method)
            return (run.time,run.bytes)
        end
        a,ap,b,bp = positiveP(Xoshiro(200),s,N;method)
        n1,n2 = ap.*a,bp.*b
        D1,D2 = mean(abs2,a.-conj.(ap)),mean(abs2,b.-conj.(bp))
        # Symmetry fixes the mean difference to zero. Restore normal ordering.
        diff2 = (n1.-n2).^2 .+ n1 .+ n2
        @printf("%.1f\t%s\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\t%.4g\t%.4g\n",
                r,method,D1,D2,standard_errors(n1)...,standard_errors(n2)...,standard_errors(a.*b)...,standard_errors(diff2)...,
                1000*median(first,timings),median(last,timings)/N)
    end
end

function benchmark_high_order(repeats=32)
    println("\nRelative complex RMSE across $repeats independent ensembles")
    println("r\tmethod\tN\tRMSE (n1)^4\tRMSE (ab)^4\tRMSE (ab)^6")
    for r in (0.5,2.0), N in (1_000,10_000,100_000), method in (:canonical,:drummond)
        s = SqueezedTwoMode(r,0)
        n,m = sinh(r)^2,sinh(r)*cosh(r)
        exact = (factorial(4)*n^4,factorial(4)*m^4,factorial(6)*m^6)
        errors = map(1:repeats) do seed
            a,ap,b,_ = positiveP(Xoshiro(1234+seed),s,N;method)
            n1,ab = a.*ap,a.*b
            estimates = (mean(x->x^4,n1),mean(x->x^4,ab),mean(x->x^6,ab))
            return estimates ./ exact .- 1
        end
        rmse = ntuple(i->sqrt(mean(e->abs2(e[i]),errors)),length(exact))
        @printf("%.1f\t%s\t%d\t%.5g\t%.5g\t%.5g\n",r,method,N,rmse...)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    benchmark_two_mode()
    benchmark_high_order()
end
