# Examples of discrete maps

export logistic, loginoise, logiwcoup, doubling, cat

scalingpetX(x::Array{Float64,1}) = x
nonlinearpetX(x::Array{Float64,1}) = x.*(1-x)
sinnonlinearpetX(x::Array{Float64,1}) = sin(2*pi*x)

# Logistic
function logisticf!(x::Array{Float64,1},a::Array{Float64,1})
  x[:] = a[:] .* x .* (1.-x)
  nothing
end
logistic(alpha::Array{Float64}) = DMap(logisticf!,(x,a)->diag(a.*(1-2x)),alpha,repmat([0. 1.],length(alpha),1))
logistic(alpha::Float64=3.8) = logistic([alpha])
#logistic(alpha::Float64) = DMap((x,a)->a*x.*(1-x),(x,a)->a*(1-2x),alpha)
#logistic(alpha::Array{Float64}) = DMap((x,a)->a.*x.*(1-x),(x,a)->diag(a.*(1-2x)),alpha,repmat([0. 1.],length(alpha),1))

# L1, L2: logistic deterministic

logisticp(alpha::F64U=3.8) = Peturbation(logistic(alpha),scalingpetX)
logistic1(alpha::F64U=3.8;largs...) = IterationSchema(logisticp(alpha),"L1",logiA;largs...)
logistic2(alpha::F64U=3.8;largs...) = IterationSchema(logisticp(alpha),"L2",logiA2;largs...)
logistich(alpha::F64U=3.8,phase::Float64=0.;largs...) = IterationSchema(logisticp(alpha),"Lh",sin100A(phase);largs...) # h for hecto

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

doubling(dim::Integer=1) = DMap(doublingf!,(x,a)->2*ones(dim),(),repmat([0. 1.],dim,1),dim,true)
#doubling() = DMap(doublingf!,(x,a)->2,(),[0. 1.],1,true)
#doubling() = doubling(1)

doublingp() = Peturbation(doubling(),nonlinearpetX)
doubling1(;largs...) = IterationSchema(doublingp(),"D1",trigA;largs...)

doublingpp() = Peturbation(doubling(),sinnonlinearpetX)
doubling2(;largs...) = IterationSchema(doublingpp(),"D2",trigA;largs...)


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
