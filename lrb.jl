@everywhere setupcode = quote
  include("DiscreteMaps.jl")
  using DiscreteMaps, HDF5, JLD, Dates
  endtime = DateTime(2015,06,30,10,00,00)

  function peturbsample(M,N,deps,phase)
    DiscreteMaps.timedsample("Lh",endtime=endtime,NP=M,NCycles=1,startstring="results/lrb/rb-$(M)-$(N)-$(deps)-$(round(phase/2pi,3))",
                             Itargs=(3.8,phase),Itkwargs=DiscreteMaps.kw(samplefn=DiscreteMaps.evensamplefn,samplefnargs=(deps),
                                                                         N=N))
  end
end

@everywhere eval(setupcode)

DiscreteMaps.newpath("results/lrb")
(length(ARGS) == 1) ? (M = int(ARGS[1])) : (M = 20)
 while (now() < endtime)
  for N in [80000], deps in ([1:100]*1e-7), phase in (2pi * [0:29]/30)
    peturbsample(M,N,deps,phase)
  end
 end

