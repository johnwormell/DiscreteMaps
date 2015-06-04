#Pkg.add("ApproxFun")
#Pkg.add("Distributions")
#Pkg.update()
module DiscreteMaps

#using ApproxFun
using Distributions
using HDF5, JLD, Dates

# General fluff
F64U = Union(Float64,Array{Float64})

tomorrowmorning() = Dates.DateTime(Dates.Date(Dates.now() + Dates.Hour(16)))+Dates.Hour(10)
newpath(path) = ispath(path) || mkdir(path)
restrictto(x::Array{Float64,1},xmin,xmax) = max(min(x,xmax),xmin)

function restrictto!(x::Array{Float64,1},xmax,xmin)
  x[:] = restrictto(x)
  nothing
end

function chopeps(epsv::Array{Float64,1},eA::Array{Float64,2},vA::Array{Float64,2};epsmax=Inf,epsmin=-Inf)
  if (epsmax < Inf )| (epsmin > -Inf)
    epsvinds = ((epsv .>= epsmin)&(epsv .<= epsmax));
    epsv = epsv[epsvinds]
    eA = eA[:,epsvinds]
    vA = vA[:,epsvinds]
  end
  return epsv, eA, vA
end

# Types

include("DM-Types.jl")

# Map iteration

include("DM-Iteration.jl")
include("DM-NewIteration.jl")

# Fluctuation-dissipation

include("DM-Fluctuation.jl")


# Noise in observables

export observe, observevar

function observe(A::Function,x_history::Array{Float64})
  # Estimates an observable
  return mean(A(x_history))::Float64#::Array{Float64} #, var(A(x_history))/length(x_history)::Array{Float64}
end

function autocovariance(A::Function,x_history::Array{Float64},sumN = 60)
  NH = length(x_history)
  A1 = A(x_history)[:]
  A2 = A(x_history)[:]
  varA = Array(Float64,sumN+1)
  for i = 0:sumN
    varA[i+1] = cov(A1,A2) #(sum(A1.*A2)-sum(A1)*sum(A2)/(N-i))/(N-i-1)
    #        flucterms[i+1] = mean(flucintegrand)
    #        flucvar[i+1] = var(flucintegrand)/(N-i)
    shift!(A1)
    pop!(A2)
  end
  return varA
end
function observevar(A::Function,x_history::Array{Float64},sumN = 60)
  # Estimates the random variance of a measured observable
  varA = autocovariance(A,x_history,sumN)
#  println(cov(A1,sort(A1)))
#  println(floor(varA))
  varA *= 2
  varA[1] /= 2
  cesaroseries = cumsum(varA)
  return sum(cesaroseries)/(sumN+1)
end

# Common discrete maps

include("DM-Examples.jl")

# Common observables

include("DM-Observables.jl")

end
