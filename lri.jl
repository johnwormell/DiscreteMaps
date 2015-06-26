println("Initialising!")
@everywhere setup_code = quote
  include("DiscreteMaps.jl")

  using DiscreteMaps, Dates

  if (length(ARGS) == 4) &&(ARGS[4] == "const")
    It = DiscreteMaps.itdict[ARGS[1]](N=int(ARGS[2]),NH=int(ARGS[3]),samplefn=DiscreteMaps.zerosamplefn)
    startstring = "results/rp"
    NQ = div(10^7,int(ARGS[2]))
  else
    startstring = "results/rs"
    NQ = 3
    if length(ARGS) == 3
      It = DiscreteMaps.itdict[ARGS[1]](N=int(ARGS[2]),NH=int(ARGS[3]))
    elseif length(ARGS) == 2
      It = DiscreteMaps.itdict[ARGS[1]](N=int(ARGS[2]))
    else
      It = DiscreteMaps.itdict[ARGS[1]]()
    end
  end

  endtime = DateTime(2015,06,29,10,00,00) #DiscreteMaps.tomorrowmorning()
end

@everywhere eval(setup_code)

println("Starting!")
DiscreteMaps.newpath("results")
cyclecount, epsv, eA, vA = DiscreteMaps.timedsample(It,endtime=endtime,NQ=NQ,startstring=startstring)

println(cyclecount)
println(epsv)
println(eA)
println(vA)
