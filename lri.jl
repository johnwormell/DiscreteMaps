println("Initialising!")
@everywhere setup_code = quote
  include("DiscreteMaps.jl")

  using DiscreteMaps, Dates

  if length(ARGS) == 3
    It = DiscreteMaps.itdict[ARGS[1]](N=int(ARGS[2]),NH=int(ARGS[3]))
  elseif length(ARGS) == 2
    It = DiscreteMaps.itdict[ARGS[1]](N=int(ARGS[2]))
  else
    It = DiscreteMaps.itdict[ARGS[1]]()
  end

  endtime = DateTime(2015,06,05,13,00,00)

end

@everywhere eval(setup_code)

println("Starting!")
cyclecount, epsv, eA, vA = DiscreteMaps.timedsample(It,endtime=endtime,NQ=3)

println(cyclecount)
println(epsv)
println(eA)
println(vA)
