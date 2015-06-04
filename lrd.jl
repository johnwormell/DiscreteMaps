@everywhere setupcode = quote
  include("DiscreteMaps.jl")
  using DiscreteMaps, HDF5, JLD, Dates
  NI = 10^5
  N = 10^5
  sumN = 20
  epsr = [0,0.0001,0.005]
  epsN = length(epsr)
  starttime = now()
  DiscreteMaps.newpath("results")
  startstring = "results/rd"
  function getdiagnostics(PI)
    It = DiscreteMaps.itdict[PI](NI=NI,N=N,NH=N)
    filename = replace("$(startstring)-$(It.PInitial)-$(starttime).h5",":","-")

    println(("Starting $PI"))
    vAv = Array(Float64,sumN,epsN)
    Av = Array(Float64,N,epsN)
    for i = 1:epsN
      xh = DiscreteMaps.obsiterate(It,epsr[i],returnhistory=true)
      vAv[:,i] = DiscreteMaps.autocovariance(It.A[1],xh,sumN)
      Av[:,i] = It.A[1](xh)
    end

    save(filename,"vAv",vAv,"Av",Av)
    println("File saved for $PI")
    return "A"
  end
end

@everywhere eval(setupcode)
println("Eval done")
pmap(getdiagnostics,ARGS)
