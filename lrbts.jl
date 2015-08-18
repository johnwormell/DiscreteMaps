@everywhere include("DiscreteMaps.jl")

@everywhere module expas
    using DiscreteMaps
    alpha = 3.8
    copts = DiscreteMaps.logisticcriticalorbit(DiscreteMaps.logistic(alpha)).pts |> vec

    function trigo(center::Float64,k::Int64)
        A(x::Array{Float64}) = cos(2*pi*k*(x[1,:]-center))
    end
    trigcoA(kv,ns) = [DiscreteMaps.ido(),vec([trigo(copts[i],k) for k in kv,i in ns])]
    logisticts(kv,ns;largs...) = DiscreteMaps.IterationSchema(DiscreteMaps.logisticp(alpha),"Lts",trigcoA(kv,ns);largs...)
    end

@everywhere setupcode = quote
    nv = [1,5,13,18,22]
    using DiscreteMaps


    using expas
    function peturbsample(M,deps,kr)
        N = 10^8
        NH = 10^2
        kv = (sqrt(3) / 2pi) ./ 2. .^(-[1:0.5:27]) |> round |> int
        DiscreteMaps.timedsample(expas.logisticts(kv,nv,
            samplefn=DiscreteMaps.evensamplefn,samplefnargs=(deps),N=N,NH=NH),
            NP=M,NCycles=1,startstring="results/lrb/rbugs-$(deps)-$(kr)-")
    end
end
@everywhere eval(setupcode)
M = 20
deps = 1e-6
peturbsample(M,deps,1.);

