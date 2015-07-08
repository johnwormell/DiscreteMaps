#Pkg.add("ApproxFun")
#Pkg.add("Distributions")
#Pkg.add("Roots")
#Pkg.update()
module DiscreteMaps

#using ApproxFun
using Distributions
using HDF5, JLD, Dates
using Roots

# General fluff

# Input types
F64U = Union(Float64,Array{Float64})
I64U = Union(Int64,Array{Int64})

# returns 10AM tomorrow
tomorrowmorning() = Dates.DateTime(Dates.Date(Dates.now() + Dates.Hour(16)))+Dates.Hour(10)

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
#include("DM-Display.jl")

## SPECTRAL APPROXIMATION OF INVARIANT MEASURES
# Spectral functions
include("DM-Spectral.jl")

# Measures
include("DM-Measures.jl")

# Spectral acim
include("DM-Acim.jl")

end
