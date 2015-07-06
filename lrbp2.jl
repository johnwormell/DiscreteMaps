@everywhere setupcode = quote
  include("DiscreteMaps.jl")
  using DiscreteMaps, HDF5, JLD, Dates
  endtime = DateTime(2016,07,06,12,00,00)
  function peturbsample(M,deps)
    N = 10^7
    NH = 10^2
    DiscreteMaps.timedsample("Lh",NP=M,NCycles=1,startstring="results/lrb/rbp-$(deps)",
                             Itargs=(3.8),Itkwargs=DiscreteMaps.kw(samplefn=DiscreteMaps.evensamplefn,samplefnargs=(deps),N=N,NH=NH))
  end
end

@everywhere eval(setupcode)

DiscreteMaps.newpath("results/")
DiscreteMaps.newpath("results/lrb")
(length(ARGS) == 1) ? (M = int(ARGS[1])) : (M = 20)
 while (now() < endtime)
  for deps in ([1:200]*2.5e-8)
    peturbsample(M,deps)
  end
end

