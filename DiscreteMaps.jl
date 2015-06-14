#Pkg.add("ApproxFun")
#Pkg.add("Distributions")
#Pkg.update()
module DiscreteMaps

#using ApproxFun
using Distributions
using HDF5, JLD, Dates

# General fluff
F64U = Union(Float64,Array{Float64})
I64U = Union(Int64,Array{Int64})

# returns 10AM tomorrow
tomorrowmorning() = Dates.DateTime(Dates.Date(Dates.now() + Dates.Hour(16)))+Dates.Hour(10)

# Create a folder if it doesn't already exist
newpath(path) = ispath(path) || mkdir(path)

# Domain stuff
restrictto(x::Array{Float64,1},xmin,xmax) = max(min(x,xmax),xmin)

function restrictto!(x::Array{Float64,1},xmax,xmin)
  x[:] = restrictto(x)
  nothing
end

function indomain(x::Array{Float64},dom::Array{Float64,2})
  return minimum(dom[:,1] .<= x .<=dom[:,2],1)
end

function checkindomain(x::Array{Float64},dom::Array{Float64,2})
  return x[:,vec(indomain(x,dom))]
end

# Distance from edge of (hyper-)rectangular domain
function domainedgedist(x::Array{Float64},dom::Array{Float64,2})
  return (minimum(min(abs(x.-dom[:,1]),abs(x.-dom[:,2])),1))
end

# Test function!
testfn(x,centre,width) = (abs(x.-centre) .< width) .* exp(1-(1 - ((x .-centre)/width).^2).^(-1))


# Types

include("DM-Types.jl")

# Observation and noise in observables

include("DM-Observing.jl")

# Map iteration

include("DM-Iteration.jl")
include("DM-NewIteration.jl")

# Fluctuation-dissipation

include("DM-Fluctuation.jl")

# Common discrete maps

include("DM-Examples.jl")

# Common observables

include("DM-Observables.jl")

# Display

include("DM-Display.jl")

# ACIM

include("DM-Acim.jl")

end
