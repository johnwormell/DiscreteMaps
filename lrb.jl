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
for i = 1:20
  peturbsample(20,40000,1e-3)
  peturbsample(20,40000,1e-4)
  peturbsample(20,40000,1e-5)
  peturbsample(20,40000,1e-6)
  peturbsample(40,40000,1e-3)
  peturbsample(40,40000,1e-4)
  peturbsample(40,40000,1e-5)
  peturbsample(40,40000,1e-6)
  peturbsample(20,80000,1e-3)
  peturbsample(20,80000,1e-4)
  peturbsample(20,80000,1e-5)
  peturbsample(20,80000,1e-6)
  peturbsample(40,80000,1e-3)
  peturbsample(40,80000,1e-4)
  peturbsample(40,80000,1e-5)
  peturbsample(40,80000,1e-6)
end

