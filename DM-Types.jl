export Map, NDMap, DMap, Peturbation

# Map
abstract Map


type NDMap <: Map
  f::Function
  params
  dom::Array
  dim::Integer
  periodic::Bool
  init::Function
  noise::Function
end

runifinit(dom,dim) = return dom[:,1]+rand(dim).*(dom[:,2]-dom[:,1])

function rnormadditivenoise(dom,dim)
  function noise!(x::Array{Float64})
    pet = 100*(dom[:,2]-dom[:,1])*eps(Float64).*randn(dim)
    x[:] = max(min(x+pet,1),0)
    nothing
  end
  return noise!
end

NDMap(f!,params,dom,dim) = NDMap(f!,params,dom,dim,false,()->runifinit(dom,dim),rnormadditivenoise(dom,dim))
NDMap(f!,params,dom) = NDMap(f!,params,dom,div(length(dom),2))
NDMap(f!,params) = NDMap(f!,params,[0. 1.])

# DifferentiableMap
type DMap <: Map
  f!::Function
  df::Function
  params
  dom::Array
  dim::Integer
  periodic::Bool
  init::Function
  noise!::Function
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
  samplefnargs
  N::Integer
  NI::Integer
  NH::Integer
  AN::Integer
end
IterationSchema(P::Peturbation,Pinitial::String,A::Array{Function,1};
                samplefn::Function=betasamplefn,samplefnargs=(),
                N::Integer=10^7,NI::Integer=10^4,NH::Integer=10^4) =
  IterationSchema(P,Pinitial,A,samplefn,samplefnargs,N,NI,NH,length(A))

