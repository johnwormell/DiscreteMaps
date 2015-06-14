export Map, NDMap, DMap, Peturbation

runifinit(dom,dim) = return dom[:,1]+rand(dim).*(dom[:,2]-dom[:,1])

function rnormadditivenoise(dom,dim)
  function noise!(x::Array{Float64})
    pet = 100*(dom[:,2]-dom[:,1])*eps(Float64).*randn(dim)
    x[:] = max(min(x+pet,1),0)
    nothing
  end
  return noise!
end

function makef(f!::Function)
  function f(x::Array{Float64,1},a)
    P = copy(x)
    f!(P,a)
    return P
  end
  f(x::Float64,a) = f([x],a)
  return f
end

# Map
abstract Map

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

# DifferentiableMap
type DMap <: Map
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

# type Acim
#    M::DMap
#    rho::Function
#    sp::Space
# end

DMap(M::Map,df) = DMap(M.f!,df,M.params,M.dom,M.dim,M.periodic,M.init,M.noise);
DMap(f!,df,params,dom,dim,periodic) = DMap(f!,df,params,dom,dim,periodic,()->runifinit(dom,dim),rnormadditivenoise(dom,dim))
DMap(f!,df,params,dom,dim) = DMap(f!,df,params,dom,dim,false,()->runifinit(dom,dim),rnormadditivenoise(dom,dim))
DMap(f!,df,params,dom) = DMap(f!,df,params,dom,div(length(dom),2))
DMap(f!,df,params) = DMap(f!,df,params,[0. 1.])
# DMap(P::Peturbation,epsilon::)
#TODO Make this ^^

type IMap <: Map
  f::Function
  f!::Function
  df::Function
  g::Function
  params
  crit::Function
  critxx::Function
  dom::Array
  dim::Integer
  periodic::Bool
  init::Function
  noise!::Function

  IMap(f!,df,g,params,crit,critxx,dom,dim,periodic,init,noise!) = new(makef(f!),f!,df,g,params,crit,critxx,dom,dim,periodic,init,noise!)
#  IMap(f!,args...) = addf(new,f!,args...)
end

IMap(M::DMap,g,crit,critxx) = IMap(M.f!,M.df,g,M.params,crit,critxx,M.dom,M.dim,M.periodic,M.init,M.noise!);
#DMap(f!,df,params,dom,dim,periodic) = DMap(f!,df,params,dom,dim,periodic,()->runifinit(dom,dim),rnormadditivenoise(dom,dim))
#DMap(f!,df,params,dom,dim) = DMap(f!,df,params,dom,dim,false,()->runifinit(dom,dim),rnormadditivenoise(dom,dim))
#DMap(f!,df,params,dom) = DMap(f!,df,params,dom,div(length(dom),2))
#DMap(f!,df,params) = DMap(f!,df,params,[0. 1.])

# Peturbation

type Peturbation
  M::Map
  X::Function
  defaulteps::Float64
end

Peturbation(M::Map,X) = Peturbation(M::Map,X,0.05)

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

type CriticalOrbit
  crit::Array{Float64,1}
  pts::Array{Float64,2}
  mag::Array{Float64,2}
  sgn::Array{Int8,2}
  Nc::Integer
  Npts::Integer
end
CriticalOrbit(crit,pts,mag,sign) = CriticalOrbit(crit,pts,mag,sign,length(crit),size(pts,1))

type Spikes
  CO::CriticalOrbit
  dom::Array{Float64,2}
  mag0::Array{Float64,1}
  width::Array{Float64,1}
  bumpwidth::Array{Float64,1}
end
Spikes(CO::CriticalOrbit,dom::Array{Float64,2},mag0=fill(0.01,length(CO.crit)),
       width::Float64=fill(0.01,length(CO.crit)),
       bumpwidth::Array{Float64,1} = fill(0.05,length(CO.crit)),
       ) =
                     Spikes(CO,dom,mag0,width,bumpwidth)

# type Acim
#    M::DMap
#    rho::Function
#    sp::Space
# end
