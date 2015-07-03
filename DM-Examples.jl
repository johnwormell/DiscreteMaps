# Examples of discrete maps

scalingpetX(x::Array{Float64,1}) = x
nonlinearpetX(x::Array{Float64,1}) = x.*(1-x)
sinnonlinearpetX(x::Array{Float64,1}) = sin(2*pi*x)

# Logistic

# So logisticg etc. won't throw a DomainError if x isn't in the image of f,
# but just output empty vector

function nicesqrt(x::F64U)
  minimum(x) >= 0 ? (sqrt(x)) : (fill(convert(typeof(x[1]),NaN),size(x)))
end

# Map functions
function logisticf!(x::Array{Float64,1},a::Array{Float64,1})
  x[:] = a[:] .* x .* (1.-x)
  nothing
end

function makelogisticg(;newdom::Bool=false)
  function logisticg(x::F64U,a::F64U)
    g1 = (1 - nicesqrt(1 - 4*x./a))/2 # + eps(maximum(a)/4)
    g2 = 1 - g1
    return checkindomain([g1 g2],logisticdom(a;newdom=newdom))
  end
#  logisticg(x::Float64,a) = logisticg([x]::Array{Float64,1},a::Array{Float64,1}) # SKULL AND XBONES FOR PARALLELISATION
  return logisticg
end

# Logistic domain
# This is not the natural domain of the logistic map in any sense but it
# limits the number of discontinuities when you do the transfer map so the
# number of artefact functions required, and for that we must be greatful
function logisticdom(alpha::Array{Float64,1};newdom::Bool=false)
  if newdom
    return [alpha.^2.*(4-alpha)/16 alpha/4]
  else
    return [1-alpha/4 alpha/4]
  end
end
logisticdom(alpha::Float64;newdom::Bool=false) = logisticdom([alpha];newdom=newdom)

# Logistic inversetransferwt - see type definition
function logisticinversetransferwt(x::Array{Float64,1},a::Array{Float64,1})
  return fill(0.5,size(x))
#   return (x .< 0.5) .* (1 - 0.5 *testfn(x,0.5,0.1)) +
#     (x .> 0.5) .* (0.5 * testfn(x,0.5,0.1)) + # in the chaotic regime, a width of 0.1 will keep it clear of the edge of the domain for one application of inversetransfer
#     (x .== 0.5) .* 1.
# # This is really here "just in case"
end

# returns logistic Artefacts container for parameter alpha
function logisticartefacts(alpha;newdom=false)
  discont = makef(logisticf!)(logisticdom(alpha,newdom=newdom)[:,1],alpha) # where is the discontinuity happening?
  logartefn(x::F64U) = 1. * (x .< [discont]') # step function, returns a matrix
  if newdom
   logpointsin= [discont - 10*eps(maxabs(discont)), discont + 10*eps(maxabs(discont))] # get values near the edge of the step function
   loggetcoeffs(y::F64U,i::Int64) = y[2] - y[1] # we use the difference in values as a coeff
  else
    logpointsin = discont + 2*eps(maxabs(discont)) # get value at the edge of the step function
    loggetcoeffs(y::F64U,i::Int64) = y # we just use the value
  end
  lognfns = 1 # at present

  return Artefacts(logartefn,logpointsin,loggetcoeffs,lognfns)
end

function logistic(alpha::Array{Float64},invertible=(length(alpha)==1);newdom::Bool=false)
  if ~(invertible)
    DMap(logisticf!, # logistic map
         (x,a) -> (a .* (1 - 2x)), # differential of 1D logistic map
         alpha, # parameter is the logistic parameter
         logisticdom(alpha;newdom=newdom)) # domain of logistic map
    # the rest is specified by defaults in DMap
  else
    IMap(logistic(alpha,false;newdom=newdom), # logistic DMap, see immediately above
         makelogisticg(newdom=newdom), # inverse of logistic map
         logisticinversetransferwt, # logisticinversetransferwt, just constant 1/2
         (a)->[0.5], # one critical point
         (a)->-2*a, # f'(0.5) = -2α
         logisticartefacts(alpha;newdom=newdom) # logistic Artefacts container
         )
  end
end
logistic(alpha::Float64=3.8;newdom=false) = logistic([alpha],true;newdom=newdom)
#logistic(alpha::Float64) = DMap((x,a)->a*x.*(1-x),(x,a)->a*(1-2x),alpha)
#logistic(alpha::Array{Float64}) = DMap((x,a)->a.*x.*(1-x),(x,a)->diagm(a.*(1-2x)),alpha,repmat([0. 1.],length(alpha),1))

### L1, L2, Lh: logistic deterministic

logisticp(alpha::F64U=3.8) = Peturbation(logistic(alpha),scalingpetX)
logistic1(alpha::F64U=3.8;largs...) = IterationSchema(logisticp(alpha),"L1",logiA;largs...)
logistic2(alpha::F64U=3.8;largs...) = IterationSchema(logisticp(alpha),"L2",logiA2;largs...)
logistich(alpha::F64U=3.8,phase::Float64=0.;largs...) = IterationSchema(logisticp(alpha),"Lh",sin100LhdefaultcptsA(phase);largs...) # h for hecto
logistichp(alpha::F64U=3.8;largs...) = IterationSchema(logisticp(alpha),"Lhp",sin100p30A;largs...) # h for hecto
# Logistic with noise
function loginoisef!(x::Array{Float64,1},a::(Array{Float64,1},Float64))
  x[:] = restrictto(a[1] .* x .* (1.-x) + a[2]*randn(length(x)),0.,1.)
  nothing
end

loginoise(alpha::Array{Float64,1},sd::Float64) = DMap(loginoisef!,(x,a)->max(min(a[1]+a[2]*randn(1),4.),0.).*(1-2x),(alpha,sd))
loginoise(alpha::Float64=3.8,sd::Float64=0.02) = loginoise([alpha],sd)

### N1: big noise

loginoisep(alpha::F64U=3.8,sd::Float64=0.02) = Peturbation(loginoise(alpha,sd),scalingpetX)
loginoise1(alpha::F64U=3.8,sd::Float64=0.02;largs...) = IterationSchema(loginoisep(alpha,sd),"N1",logiA;largs...)

### M1: medium noise

loginoisem1(alpha::F64U=3.8,sd::Float64=0.0002;largs...) = IterationSchema(loginoisep(alpha,sd),"M1",logiA;largs...)

# Doubling

function doublingf!(x::Array{Float64,1},a::())
  x[:] = mod(2*x,1)
  nothing
end

function doublingg(x::Array{Float64,1},a::())
  g1 = x/2
  g2 = (x+1)/2
  return [g1 g2]
end
doublingg(x::Float64,a) = doublingg([x],a)

#doubling(dim::Integer=1) = DMap(doublingf!,(x,a)->2*ones(size(x)),(),repmat([0. 1.],dim,1),dim,true)
#doubling() = DMap(doublingf!,(x,a)->2,(),[0. 1.],1,true)
#doubling() = doubling(1)

function doubling(dim::Integer=1,invertible::Bool=(dim==1))
  if ~(invertible)
        DMap(doublingf!, # double it!
             (x,a)->2*ones(size(x)), # let its derivative be 2!
             (), # no parameters
             repmat([0. 1.],dim,1), # got a unit domain!
             dim, # of an arbitrary number of dimensions but by default 1 (see above)
             true) # it's periodic
  else
    IMap(doubling(dim,false), #doubling DMap, see immediately above
         doublingg, # inverse of doubling map
         (a)->[],(a)->[], # critical points stuff (none)
         (x,a) -> fill(0.5,size(x)), # inversetransferwt fn - evenly distributed
         Artefacts()) # no artefacts

  end
end


doublingp() = Peturbation(doubling(),nonlinearpetX)
doubling1(;largs...) = IterationSchema(doublingp(),"D1",trigA;largs...)

doublingpp() = Peturbation(doubling(),sinnonlinearpetX)
doubling2(;largs...) = IterationSchema(doublingpp(),"D2",trigA;largs...)

# Sine-peturbed doubling

spdoublinglift(x::Float64,a::Float64) = 2*(x + a*sin(2pi*x))
function spdoublingf!(x::Array{Float64,1},a::Float64)
  x[:] = mod(spdoublinglift(x[1],a),1)
  nothing
end

function spdoublingg(y::Array{Float64,1},a::Float64)
  g1 = fzero(x->spdoublinglift(x,a) - y[1],0.,1.)
  g2 = fzero(x->spdoublinglift(x,a) - (y[1]+1),0.,1.)
  return [g1 g2]
end
spdoublingg(x::Float64,a) = spdoublingg([x],a)

function spdoubling(a::Float64=0.05,invertible::Bool=true)
  if ~(invertible)
        DMap(spdoublingf!, # double it + peturbation
             (x,a)->2*(1 + 2pi*a*cos(2pi*x)), # derivative
             a, # a peturbation parameter
             repmat([0. 1.],1,1), # got a unit domain of dimension 1
             1, # one dimension remember
             true) # it's periodic
  else
    IMap(spdoubling(a,false), #spdoubling DMap, see immediately above
         spdoublingg, # inverse of spdoubling map
         (x,a) -> fill(0.5,size(x)) # inversetransferwt fn - evenly distributed
         ) # no artefacts or critical points
  end
end


# Arnol'd cat map

catmatrix = [2. 1.;1. 1.]
function catf!(x::Array{Float64,1},a::())
  x[:] = mod(catmatrix*x,1)
end

catm() = DMap(catf!,(x,a)->catmatrix,(),[0. 1.;0. 1.],2,true)
catp() = Peturbation(catm(),nonlinearpetX)
cat1(;largs...) = IterationSchema(catp(),"C1",trig2A;largs...)


catpp() = Peturbation(catm(),sinnonlinearpetX)
cat2(;largs...) = IterationSchema(catpp(),"C2",trig2A;largs...)

# Logistic coupling

logicoup_choosealphas(J::Integer,alpha0::Float64,sd0::Float64) = restrictto(alpha0+sd0*randn(J),0,4)

## no coupling
loginocoupp(J::Integer=100,alpha0::Float64=3.8,sd0::Float64=0.02) =
  logisticp(logicoup_choosealphas(J,alpha0,sd0))
loginocoup1(J::Integer=100,alpha0::Float64=3.8,sd0::Float64=0.02;largs...) =
  IterationSchema(loginocoupp(J,alpha0,sd0),"Y1",coupA;largs...)

## with coupling
function logiwcmap!(x::Array{Float64,1},a::(Array{Float64,1},Float64,Int64,Int64))
  alpha = a[1]
  f, Jind = gj(a[3],a[4])
  Jmult = length(alpha)^(Jind)
  phit = Jmult * f(x)[1]
  x[:] = (alpha + a[2] * phit) .* x .* (1-x)
end

function logiwcoup(alpha::Array{Float64},sd::Float64,n::Int64,gj::Int64)
  a = (alpha,sd,n,gj)
  return DMap(logiwcmap!,identity, # fix
              a, repmat([0. 1.],length(alpha),1))
end

### W1: georg

logiwcoupgp(J::Integer=100;alpha0::Float64=3.8,sd0::Float64=0.02,n::Integer=100,sd::Float64=0.02) =
  Peturbation(logiwcoup(logicoup_choosealphas(J,alpha0,sd0),sd,n,1),scalingpetX) # fix
logiwcoupg1(J::Integer=100;alpha0::Float64=3.8,sd0::Float64=0.02,n::Integer=100,sd::Float64=0.02,largs...) =
  IterationSchema(logiwcoupgp(J,alpha0=alpha0,sd0=sd0,n=n,sd=sd),"W1",coupA;largs...) # fix

### W2: gewgaw

logiwcoupggp(J::Integer=100;alpha0::Float64=3.8,sd0::Float64=0.02,n::Integer=100,sd::Float64=0.02) =
  Peturbation(logiwcoup(logicoup_choosealphas(J,alpha0,sd0),sd,n,3),scalingpetX) # fix
logiwcoupg2(J::Integer=100;alpha0::Float64=3.8,sd0::Float64=0.02,n::Integer=100,sd::Float64=0.02,largs...) =
  IterationSchema(logiwcoupggp(J,alpha0=alpha0,sd0=sd0,n=n,sd=sd),"W2",coupA2;largs...) # fix

### W3: gewgaw + different scaling/gewgawaux

logiwcoupgggp(J::Integer=100;alpha0::Float64=3.8,sd0::Float64=0.02,n::Integer=100,sd::Float64=0.02) =
  Peturbation(logiwcoup(logicoup_choosealphas(J,alpha0,sd0),sd,n,4),scalingpetX) # fix
logiwcoupg3(J::Integer=100;alpha0::Float64=3.8,sd0::Float64=0.02,n::Integer=100,sd::Float64=0.02,largs...) =
  IterationSchema(logiwcoupgggp(J,alpha0=alpha0,sd0=sd0,n=n,sd=sd),"W3",coupA2;largs...) # fix


### X1: jeroen
logiwcoupjp(J::Integer=100;alpha0::Float64=3.8,sd0::Float64=0.02,n::Integer=100,sd::Float64=0.02) =
  Peturbation(logiwcoup(logicoup_choosealphas(J,alpha0,sd0),sd,n,2),scalingpetX) # fix
logiwcoupj1(J::Integer=100;alpha0::Float64=3.8,sd0::Float64=0.02,n::Integer=100,sd::Float64=0.02,largs...) =
  IterationSchema(logiwcoupjp(J,alpha0=alpha0,sd0=sd0,n=n,sd=sd),"X1",coupA;largs...) # fix

### X2: jinx

logiwcoupjnp(J::Integer=100;alpha0::Float64=3.8,sd0::Float64=0.02,n::Integer=100,sd::Float64=0.02) =
  Peturbation(logiwcoup(logicoup_choosealphas(J,alpha0,sd0),sd,n,4),scalingpetX) # fix
logiwcoupj2(J::Integer=100;alpha0::Float64=3.8,sd0::Float64=0.02,n::Integer=100,sd::Float64=0.02,largs...) =
  IterationSchema(logiwcoupjnp(J,alpha0=alpha0,sd0=sd0,n=n,sd=sd),"X2",coupA2;largs...) # fix

### X3: jinx + different scaling

logiwcoupjnnp(J::Integer=100;alpha0::Float64=3.8,sd0::Float64=0.02,n::Integer=100,sd::Float64=0.02) =
  Peturbation(logiwcoup(logicoup_choosealphas(J,alpha0,sd0),sd,n,6),scalingpetX) # fix
logiwcoupj3(J::Integer=100;alpha0::Float64=3.8,sd0::Float64=0.02,n::Integer=100,sd::Float64=0.02,largs...) =
  IterationSchema(logiwcoupjnnp(J,alpha0=alpha0,sd0=sd0,n=n,sd=sd),"X3",coupA2;largs...) # fix


# Dictionary

itdict = {"L1" => logistic1,
          "L2" => logistic2,
          "Lh" => logistich,
          "Lhp" => logistichp,
          "M1" => loginoisem1,
          "N1" => loginoise1,
          "D1" => doubling1,
          "D2" => doubling2,
          "C1" => cat1,
          "C2" => cat2,
          "Y1" => loginocoup1,
          "W1" => logiwcoupg1,
          "W2" => logiwcoupg2,
          "W3" => logiwcoupg3,
          "X1" => logiwcoupj1,
          "X2" => logiwcoupj2,
          "X3" => logiwcoupj3
          }
