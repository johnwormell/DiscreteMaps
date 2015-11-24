module DiscreteMaps

#using ApproxFun
using Distributions

using Roots

# General fluff

# Input types
F64U = Union(Float64,Array{Float64})
I64U = Union(Int64,Array{Int64})

# turn into keyword argument dictionary
kw(;kwargs...) = kwargs

# returns 10AM tomorrow
tomorrowmorning() =
  Dates.DateTime(Dates.Date(Dates.now() + Dates.Hour(16))) +
  Dates.Hour(10)

# Create a folder if it doesn't already exist
newpath(path) = ispath(path) || mkdir(path)

# Naturally extending ones
Base.ones(n::()) = 1.


## GENERAL
# Types
include("DM-Types.jl")

# Domain
include("DM-Domain.jl")

# Common discrete maps
include("DM-Examples.jl")

# Common observables
include("DM-Observables.jl")

## ERGODIC APPROXIMATION OF INVARIANT MEASURES
# Observation and noise in observables
include("DM-Observing.jl")

# Map iteration
include("DM-Iteration.jl")
include("DM-NewIteration.jl")

# Fluctuation-dissipation
include("DM-Fluctuation.jl")

# Displaying iteration stuff
include("DM-Display.jl")

## SPECTRAL APPROXIMATION OF INVARIANT MEASURES
# Spectral functions
include("DM-Spectral.jl")

# Measures
include("DM-Measures.jl")

# Spectral acim
include("DM-Acim.jl")

function logisticgaussevenAc(k::Float64;
                             granularity=0.2,alpha::Float64=3.8)
  dom = vec(logistic(alpha).dom)
  ctrs = linspace(dom...,
                  1+int(ceil((dom[2] - dom[1])/(granularity*k))))
  return [DiscreteMaps.gaussian(ctr,k) for ctr in ctrs],ctrs
end
logisticgaussevenA(args...;kwargs...) =
  logisticgaussevenAc(args...;kwargs...)[1]

function logisticgaussnearcoAc(k::Float64;
                               granularity=0.2,critpeakwdth=0.1,
                               Npts=25,alpha::Float64=3.8)
  CO = DiscreteMaps.criticalorbit(DiscreteMaps.logistic(3.8),Npts);
  pts = CO.pts
  wdths = (CO.mag[:]/critpeakwdth).^2
  ctrs = Float64[]
  for i = 1:Npts
    dom = (pts[i]-wdths[i],pts[i]+wdths[i])
    append!(ctrs,linspace(dom...,1+int(ceil(2*wdths[i]/
                                              (granularity*k)))))
  end
  return [DiscreteMaps.gaussian(ctr,k) for ctr in ctrs], ctrs
end
logisticgaussnearcoA(args...;kwargs...) =
  logisticgaussnearcoAc(args...;kwargs...)[1]

end

