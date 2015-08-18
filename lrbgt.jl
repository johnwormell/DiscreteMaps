@everywhere include("DiscreteMaps.jl")

@everywhere module expas
    using DiscreteMaps
    alpha = 3.8

    function trigo(center::Float64,k::Int64)
        A(x::Array{Float64}) = cos(2*pi*k*(x[1,:]-center))
    end

    copts = DiscreteMaps.logisticcriticalorbit(DiscreteMaps.logistic(alpha)).pts |> vec
    gtA(stdv,kv,ns) = [DiscreteMaps.ido(),vec([DiscreteMaps.gaussian(copts[i],std) for std in stdv,i in ns]),
        vec([trigo(copts[i],k) for k in kv,i in ns])]

    logisticgt(stdv,kv,ns;largs...) = DiscreteMaps.IterationSchema(DiscreteMaps.logisticp(alpha),"Lgt",gtA(stdv,kv,ns);largs...)
    end

@everywhere setupcode = quote
    nv = [1,5,13,18,22]
    using DiscreteMaps

    using expas
    function peturbsample(M,deps)
        N = 10^8
        NH = 10^2
        stdv = 2. .^(-[1:4:27])
        kv = (sqrt(3) / 2pi) ./ stdv |> round |> int
        DiscreteMaps.timedsample(expas.logisticgt(stdv,kv,nv,
            samplefn=DiscreteMaps.evensamplefn,samplefnargs=(deps),N=N,NH=NH),
            NP=M,NCycles=1,startstring="results/lrb/rbugs-$(deps)-")
    end
end
@everywhere eval(setupcode)
M = 20
for deps = [1:100]/100*1e-6
    peturbsample(M,deps);
end
