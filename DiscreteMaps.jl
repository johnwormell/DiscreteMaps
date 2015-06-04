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

# Types

include("DM-Types.jl")

# Map iteration

include("DM-Iteration.jl")
include("DM-NewIteration.jl")

# Fluctuation-dissipation

include("DM-Fluctuation.jl")

# Observation and noise in observables

include("DM-Observing.jl")

# Common discrete maps

include("DM-Examples.jl")

# Common observables

include("DM-Observables.jl")

# Display

include("DM-Display.jl")

end
