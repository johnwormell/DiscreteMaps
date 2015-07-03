@everywhere setupcode = quote
  include("DiscreteMaps.jl")
  println("include done")
  using DiscreteMaps
  using HDF5, JLD, Dates
  println("using done")
  endtime = DateTime(2015,07,06,12,00,00)
#   Nk = 10
#   Ntht = 10
#   Npts = 40

#   peakedgeheight = 0.1
#   CO = DiscreteMaps.criticalorbit(DiscreteMaps.logistic(3.8),Npts);
#   spds = DiscreteMaps.logisticcospeeds(CO,DiscreteMaps.logistic(3.8)) |> vec;
#   pts = CO.pts[:]
#   wdths = (CO.mag[:]/peakedgeheight).^2


#   relkvals = Array(Float64,1,1,Nk)
#   relkvals[:] = [1:Nk]/4
#   kvals = relkvals .* ones(1,Ntht) .* spds * 1e-5

#   relthtvals = [1:Ntht]'
#   thtvals = relthtvals .* kvals

#   sdvs = kvals[:]
#   ctrs = (pts .+ thtvals)[:]
#   typeof(sdvs[1]) |> println
#   typeof(ctrs[1]) |> println
#   NA = length(sdvs)
# #  typeof(DiscreteMaps.gaussian(ctrs[i],sdvs[i])) |> println
#   gaussA = [DiscreteMaps.gaussian(ctrs[i],sdvs[i]) for i = 1:NA]

#   for i = 1:length(ctrs)
#     A = eye
#  #   A = DiscreteMaps.gaussian(ctrs[i],sdvs[i])
#     push!(gaussA, A)
#   end

  N = 4*10^7
  NH = 4*10^4

  function peturbsample(M,deps)
    It = DiscreteMaps.IterationSchema(DiscreteMaps.logisticp(3.8),"Lg",DiscreteMaps.logisticgaussA;
                                      samplefn=DiscreteMaps.evensamplefn,samplefnargs=(deps),N=N,NH=NH)
    DiscreteMaps.timedsample(It,NP=M,NCycles=1,startstring="results/lrb/rbg-$(deps)")
  end
end
println("setup code defined")
@everywhere eval(setupcode)
print("eval done")

DiscreteMaps.newpath("results")
DiscreteMaps.newpath("results/lrb")
(length(ARGS) == 1) ? (M = int(ARGS[1])) : (M = 20)

while (now() < endtime)
  for deps in ([1:1000]*1e-8)
    peturbsample(M,deps)
  end
end

Pkg.update()

