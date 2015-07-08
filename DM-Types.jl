export Map, NDMap, DMap, Peturbation

# Domain

# type Domain
#   dom::Array{Float64,2}
#   periodic::Bool
#   dim::Integer
#   Domain(dom::Array{Float64,2},periodic::Bool) = new(dom,periodic,size(dom,1))
# end
# UnitInterval(dim=1) = Domain(repmat([0. 1.],dim,1),false)
# UnitPInterval(dim=1) = Domain(repmat([0. 1.],dim,1),true)

# Map
abstract Map
abstract DifferentiableMap <: Map

# Initial random value
runifinit(dom::Array{Float64,2},dim::Integer) = dom[:,1]+rand(dim).*(dom[:,2]-dom[:,1])
#runifinit(D:Domain) = runifinit(D.dom,D.dim)

# Standard noise function to stop iteration being periodic
function rnormadditivenoise(dom::Array{Float64,2},dim::Integer)
  function noise!(x::Array{Float64})
    pet = 100*(dom[:,2]-dom[:,1])*eps(Float64).*randn(dim)
    x[:] = max(min(x+pet,1),0)
    nothing
  end
  return noise!
end
#rnormadditivenoise(D::Domain) = rnormadditivenoise(D.dom,D.dim)

# Turns map function that changes the input into a function that returns a value
function makef(f!::Function)
  function f(x::Array{Float64,1},a)
    P = copy(x)
    f!(P,a)
    return P
  end
  f(x::Float64,a) = f([x],a)
  return f
end

# Non-differentiable map
# For explanations see below
type NDMap <: Map
  f::Function
  f!::Function
  params
  dom::Array
  dim::Integer
  periodic::Bool
  init::Function
  noise!::Function

  NDMap(f!,params,dom,dim,periodic,init,noise!) = new(makef(f!),f!,params,dom,dim,periodic,init,noise!)
end

NDMap(f!,params,dom,dim) = NDMap(f!,params,dom,dim,false,()->runifinit(dom,dim),rnormadditivenoise(dom,dim))
NDMap(f!,params,dom) = NDMap(f!,params,dom,div(length(dom),2))
NDMap(f!,params) = NDMap(f!,params,[0. 1.])

# Differentiable map
# For explanations see below
type DMap <: DifferentiableMap
  f::Function
  f!::Function
  df::Function
  params
  dom::Array
  dim::Integer
  periodic::Bool
  init::Function
  noise!::Function

  DMap(f!,df,params,dom,dim,periodic,init,noise!) = new(makef(f!),f!,df,params,dom,dim,periodic,init,noise!)
end

DMap(M::Map,df) = DMap(M.f!,df,M.params,M.dom,M.dim,M.periodic,M.init,M.noise);
DMap(f!,df,params,dom,dim,periodic) = DMap(f!,df,params,dom,dim,periodic,()->runifinit(dom,dim),rnormadditivenoise(dom,dim))
DMap(f!,df,params,dom,dim) = DMap(f!,df,params,dom,dim,false,()->runifinit(dom,dim),rnormadditivenoise(dom,dim))
DMap(f!,df,params,dom) = DMap(f!,df,params,dom,div(length(dom),2))
DMap(f!,df,params) = DMap(f!,df,params,[0. 1.])
# DMap(P::Peturbation,epsilon::)
#TODO Make this ^^

# type Acim
#    M::DMap
#    rho::Function
#    sp::Space
# end

# Container that describes artefact functions for acim determination
type Artefacts
  artefn::Function
  # takes column vector of x-values as input,
  # outputs length(x) by nfns matrix artefact function(s) evaluated at the points
  pointsin::Array{Float64,1}
  # x values to input into transfer map: the output goes to getcoeffs
  getcoeffs::Function
  # takes L(<basis function>)(pointsout) and index of basisfn as input, outputs coefficients of artefn
  nfns::Integer
  # Number of artefact functions
end
# for when there are no artefacts here, an easy way to tell Indiana goodbye
Artefacts() = Artefacts(x->Array(Float64,length(x),0),[],x->[],0)

# Invertible, differentiable map
type IMap <: DifferentiableMap
  f::Function # the map
  f!::Function # the map, but modifies the input
  df::Function # the differential of the map. Currently only works in 1D, allows array stuff
  g::Function # the inverse of the map. Only takes 1 value as input
  params # parameters of the map
  inversetransferwt::Function # determines which right-inverse of the transfer map you use
  crit::Function # returns critical points of the map
  critxx::Function # returns f'' at critical points
  Art::Artefacts # artefact functions
  dom::Array # Domain of map
  dim::Integer # Dimension of map
  periodic::Bool # Is the map periodic??
  init::Function # Initial random noise
  noise!::Function # Noise put in at every step to stop periodicity

  IMap(f!,df,g,params,inversetransferwt,crit,critxx,Art,dom,dim,periodic,init,noise!) =
    new(makef(f!),f!,df,g,params,inversetransferwt,crit,critxx,Art,dom,dim,periodic,init,noise!) #ew
#  IMap(f!,args...) = addf(new,f!,args...)
end

IMap(M::DMap,g,inversetransferwt,crit=x->[],critxx=x->[],Art=Artefacts()) = # defaults are no critical pts, artefacts
  IMap(M.f!,M.df,g,M.params,inversetransferwt,crit,critxx,Art,M.dom,M.dim,M.periodic,M.init,M.noise!); #gross

# Peturbation
# One-parameter family of maps Mε parametrised by ε.
# Domain inherited from M
# f0(x) = M.f(x)
# fε(x) = X(M.f(x))
# Involved in linear response

type Peturbation
  M::Map # Map off which the peturbation is based
  X::Function # peturbation function: see above
  defaulteps::Float64 # default epsilon to iterate on for the peturbation
end

Peturbation(M::Map,X) = Peturbation(M::Map,X,0.0)

# IterationSchema
# Peturbation, set of observables, sampling parameters

type IterationSchema
  P::Peturbation
  PInitial::String
  A::Array{Function,1}
  samplefn::Function
  N::Integer
  NI::Integer
  NH::Integer
  AN::Integer
end
IterationSchema(P::Peturbation,Pinitial::String,A::Array{Function,1};
          samplefn::Function=betasamplefn,N::Integer=10^7,NI::Integer=10^4,NH::Integer=10^4) =
  IterationSchema(P,Pinitial,A,samplefn,N,NI,NH,length(A))

# ACIM stuff

# Critical Orbit
# Contains the forward orbit of a (set of) simple critical point(s)
# combined with a sign vector which tells you the side the forward orbit
# of points near the c.p.s go (and therefore of the spikes associated with the cps),
# and a magnitude vector which tells you the relative size of the spikes compared
# to the density at the critical point.

type CriticalOrbit
  crit::Array{Float64,1}
  pts::Array{Float64,2}
  mag::Array{Float64,2}
  sgn::Array{Int8,2}
  Nc::Integer
  Npts::Integer
end
CriticalOrbit(crit,pts,mag,sign) = CriticalOrbit(crit,pts,mag,sign,length(crit),size(pts,1))

# change the widths to keep them away from domain edges, singularities
function minwidths(CO,dom,widths,cpbd)
  widths2 = Array(typeof(widths).parameters[1],size(widths))
  for i = 1:CO.Nc
    widths2[:,i] = min(
            barrierdistoriented(CO.pts[:,i]',CO.sgn[:,i]',[dom # CO.crit
                                                 ])[:]*cpbd,
            widths[:,i])
  end
  return widths2
end

# Measures
abstract Measure

# Spikes
# Contains a set of critical orbits and details about the size of the spikes that go with them.
# Mag0 is a vector of the densities at the critical points
# Widths give the maximum supports of each family of spike
type Spikes <: Measure
  CO::CriticalOrbit
  dom::Array{Float64,2}
  mag0::Array{Float64,1}
  widths::Array{Float64,2}

#  bumpwidth::Array{Float64,1}
  Spikes(CO::CriticalOrbit,
         dom::Array{Float64,2},
         mag0::Array{Float64,1}=fill(1.,CO.Nc),
         widths::Array{Float64,2}=fill(0.05,CO.Npts,CO.Nc),
         cpbd = 0.5) =
    new(CO,dom,mag0,minwidths(CO,dom,widths,cpbd))
 #dom::Array{Float64,2}
end
(*)(normfactor::Real,Sp::Spikes) = Spikes(Sp.CO # don't change the magnitudes in this cause they're relative to mag0
                                             ,Sp.dom,normfactor*Sp.mag0,Sp.widths)

type SpectralMeasure <: Measure
  coeffs::Array{Float64,1}
  dom::Array{Float64,2}
  periodic::Bool
  N::Integer
  SpectralMeasure(coeffs::Array{Float64,1},dom::Array{Float64,2},periodic::Bool) =
    new(coeffs,dom,periodic,length(coeffs))
end
(*)(normfactor::Real,mu::SpectralMeasure) = SpectralMeasure(mu.coeffs * normfactor, mu.dom, mu.periodic)

type SumMeasure <: Measure
  components::Array{Measure}
  dom::Array{Float64,2}
  function SumMeasure(components::Array{Measure}) # can be fixed by making domain part of measure type
    dom1 = (components[1]).dom
    L = length(components)
    if length(components) > 1
      for i in [2:length(components)]
        if (components[i]).dom != dom1
          error("Measures do not have matching domain")
        end
      end
    end
    return new(components,dom1)
  end
end
function (*)(normfactor::Real,mu::SumMeasure)
  components = Measure[]
  for i = 1:length(mu.components)
    push!(components,normfactor * mu.components[i])
  end
  return SumMeasure(components)
end

(+)(mu1::Measure,mu2::Measure) = SumMeasure(Measure[mu1,mu2])
(-)(mu1::Measure,mu2::Measure) = SumMeasure(Measure[mu1,-1*mu2])
(*)(mu::Measure,normfactor::Real) = normfactor * mu
(/)(mu::Measure,normfactor::Real) = mu * (1./normfactor)
