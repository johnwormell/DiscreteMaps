if isempty(ARGS) | ~(ARGS[1] in ["cp","even"])
  error("Gaussian type must be either cp or even")
elseif ARGS[1] == "cp"
  gtype = "c"
else
  gtype = "e"
end

@everywhere setupcode = quote
  include("DiscreteMaps.jl")
  using DiscreteMaps, HDF5, JLD, Dates
  endtime = DateTime(2016,07,06,12,00,00)
  function peturbsample(M,deps,kr)
    N = 10^7
    NH = 10^2
    DiscreteMaps.timedsample("Lg$gtype",NP=M,NCycles=1,startstring="results/lrb/rbug$(gtype)-$(deps)-$(kr)-",
                             Itargs=(1/kr),Itkwargs=DiscreteMaps.kw(samplefn=DiscreteMaps.evensamplefn,samplefnargs=(deps),N=N,NH=NH))
  end
end

@everywhere eval(setupcode)

DiscreteMaps.newpath("results/")
DiscreteMaps.newpath("results/lrb")
(length(ARGS) == 2) ? (M = int(ARGS[2])) : (M = 20)
kr = 1.
krstep = 1.
while (now() < endtime)
  for deps in 10.^(linspace(-8,-5,200))
    peturbsample(M,deps,kr)
  end
  kr += krstep
end

