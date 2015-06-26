@everywhere setupcode = quote
  include("DiscreteMaps.jl")
  using DiscreteMaps, HDF5, JLD, Dates

  function peturbsample(M,N,deps)
    DiscreteMaps.timedsample("Lh",NP=M,NCycles=1,startstring="results/lrb/rb-$(M)-$(N)-$(deps)",
                             Itargs=(),Itkwargs=DiscreteMaps.kw(samplefn=DiscreteMaps.evensamplefn,samplefnargs=(deps),
                                                 N=N))
  end
end

@everywhere eval(setupcode)

DiscreteMaps.newpath("results/lrb")
while true
  for n in [40000,60000], deps in [3e-6,5e-6,1e-5,2e-5,3e-5]
    peturbsample(20,n,deps)
  end
end

