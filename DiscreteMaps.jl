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

# turn into keyword argument dictionary
kw(;kwargs...) = kwargs

# returns 10AM tomorrow
tomorrowmorning() = Dates.DateTime(Dates.Date(Dates.now() + Dates.Hour(16)))+Dates.Hour(10)

# Create a folder if it doesn't already exist
newpath(path) = ispath(path) || mkdir(path)

# Extending ones
Base.ones(n::()) = 1.

# Domain stuff
domsize(dom::Array{Float64,2}) = dom[:,2] - dom[:,1]
#domsize(M::Map) = domsize(M.dom)

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
  return broadcast(min,abs(x.-dom[:,1]),abs(x.-dom[:,2]))
end

# Distance from edge of (hyper-)rectangular domain for x pointing in oriented directions
function domainedgedistoriented(x::Array{Float64},sgn::Array{Int8},dom::Array{Float64,2})
  dommatrix = (sgn .== -1) .* dom[:,1] + (sgn .== 1) .* dom[:,2]
  return abs(x .- dommatrix)
end

# Distance from nearest (hyper-rectangular) "barrier" for x pointing in oriented directions
function barrierdistoriented(x::Array{Float64},sgn::Array{Int8},dom::Array{Float64})
  distmatrix = fill(convert(Float64,-Inf),size(x))
  pmatrix = Array(Float64,length(x))
  for i = 1:size(dom,2)
    pmatrix[:] = (dom[:,i] .- x)[:].*sgn[:]
    distmatrix[:] = max(distmatrix[:],pmatrix)
  end
  return distmatrix
end

# Test function!
testfn(x::F64U,centre,width) = (abs(x.-centre) .< width) .* exp(1-(1 - ((x .-centre)./width).^2).^(-1))
normalisedtestfnspike(x::F64U) = (x .> 0) .* testfn(x,0,1) ./ sqrt(x)
normalisedtestfnspiketotalintegral = 1.526429236397948867946243218050187830844867513079036769761039753444952626979774
  # integral of testfn(x,0,1)/sqrt(x) on [0,1]
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

include("DM-Acim2.jl")
include("DM-Acim.jl")

# logistic 3.8 gaussians - move elsewhere soon

  Nk = 6
  Ntht = 6
  Npts = 30

  peakedgeheight = 0.1
  CO = DiscreteMaps.criticalorbit(DiscreteMaps.logistic(3.8),Npts);
  spds = DiscreteMaps.logisticcospeeds(CO,DiscreteMaps.logistic(3.8)) |> vec;
  pts = CO.pts[:]
  wdths = (CO.mag[:]/peakedgeheight).^2


  relkvals = Array(Float64,1,1,Nk)
  relkvals[:] = [1:Nk]/4
  kvals = relkvals .* ones(1,Ntht) .* spds * 1e-5

  relthtvals = [1:Ntht]'
  thtvals = relthtvals .* kvals

  sdvs = kvals[:]
  ctrs = (pts .+ thtvals)[:]
  NA = length(sdvs)
#  typeof(DiscreteMaps.gaussian(ctrs[i],sdvs[i])) |> println
  logisticgaussA = [DiscreteMaps.gaussian(ctrs[i],sdvs[i]) for i = 1:NA]


end

