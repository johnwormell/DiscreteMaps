@everywhere include("DiscreteMaps.jl")

@everywhere module expas
    using DiscreteMaps
    alpha = 3.8
    copts = DiscreteMaps.logisticcriticalorbit(DiscreteMaps.logistic(alpha)).pts |> vec
    gauA(stdv,ns) = [DiscreteMaps.gaussian(copts[i],std) for std in stdv,i in ns] |> vec
    logisticgs(stdv,ns;largs...) = DiscreteMaps.IterationSchema(DiscreteMaps.logisticp(alpha),"Lgs",gauA(stdv,ns);largs...)
    end

@everywhere setupcode = quote
    nv = [1,13,22]
    using DiscreteMaps
    

    using expas
    function peturbsample(M,deps,kr)
        N = 10^8
        NH = 10^2
        DiscreteMaps.timedsample(expas.logisticgs(2. .^(-[1:0.5:30]),nv,
            samplefn=DiscreteMaps.evensamplefn,samplefnargs=(deps),N=N,NH=NH),
            NP=M,NCycles=1,startstring="results/lrb/rbugs-$(deps)-$(kr)-")
    end
end
@everywhere eval(setupcode)
peturbsample(M,deps,1.);