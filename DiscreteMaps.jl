module DiscreteMaps

using Distributions
using JLD, Dates
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
Base.zeros(n::()) = 0.

function (*){T<:Number}(t::Tridiagonal{T},x::Array{T,2})
  a = Array(Float64,size(t,1),size(x,2))
  for i = 1:size(x,2)
    a[:,i] = t*x[:,i]
  end
  a
end
# chop
chop{T<:Real}(x::T,epsl=1000eps()) =
  (abs(x) < epsl) ? convert(T,0) : x
function chopm!{T<:Real}(x::Array{T})
  epsl = 10*maximum(size(x))*eps(maxabs(x))
  for i = 1:length(x)
    x[i] = chop(x[i],epsl)
  end
  nothing
end
function chopm{T<:Real}(x::Array{T})
  epsl = 10*maximum(size(x))*eps(maxabs(x))
  xr = Array(Float64,size(x))
  for i = 1:length(x)
    xr[i] = chop(x[i],epsl)
  end
  xr
end
function chopm!{T<:Real}(x::Tridiagonal{T})
  epsl = 0.
  fields = names(Tridiagonal)
  for fld in fields
    epsl = max(epsl,100*maximum(size(x))*
                 eps(maxabs(getfield(x,fld))))
  end
  for fld in fields
    for i = 1:length(getfield(x,fld))
      getfield(x,fld)[i] = chop(getfield(x,fld)[i],epsl)
    end
  end
  nothing
end
function chopm{T<:Real}(x::Tridiagonal{T})
  xr = deepcopy(x)
  chopm!(xr)
  xr
end

## GENERAL
# Types
include("DM-Types.jl")

# Domain
include("DM-Domain.jl")

# Common discrete maps
include("DM-Examples.jl")

# Common observables
#include("DM-Observables.jl")

## ERGODIC APPROXIMATION OF INVARIANT MEASURES
# Observation and noise in observables
#include("DM-Observing.jl")

# Map iteration
include("DM-Iteration.jl")
#include("DM-NewIteration.jl")

# Fluctuation-dissipation
#include("DM-Fluctuation.jl")

# Displaying iteration stuff
include("DM-Display.jl")

## SPECTRAL APPROXIMATION OF INVARIANT MEASURES
# Spectral functions
#include("DM-Spectral.jl")
include("DM-LegSpectral.jl")

# Spectral function/spike interaction
include("DM-LegSpikes.jl")

# Measures
include("DM-Measures.jl")

# Spectral acim
include("DM-Acim.jl")

# Spectral acim function using Koopman operator rather than transfer operator
include("DM-Acimpf.jl")

end
