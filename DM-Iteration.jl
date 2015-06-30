# Iteration
export iterate, observepeturbation, checklinearresponse

function iterationpreloop!(x::Array{Float64},fn!::Function,NI::Integer)
    for n = 1:NI
        fn!(x)
    end
    nothing
end

function iterationloop!(x_history::Array{Float64},x::Array{Float64},fn!::Function,N::Integer)
    for n = 1:N
#        for j = 1:16
       fn!(x)
#        end
        x_history[:,n] = x
    end
    x
end

function iterate(fn!,M::Map,N::Integer=10^7,NI::Integer=10^4,x0::Array{Float64,1}=Float64[])
##    fn = x->DM.noise(DM.f(x,DM.params))
    x_history = Array(Float64,M.dim,N)

    if isempty(x0)
        x = M.init()
        iterationpreloop!(x,fn!,NI)
    else
    x = x0
    x0 = []
    end
    iterationloop!(x_history,x,fn!,N)
    return x_history
end

function iterate(M::Map,N::Integer=10^7,NI::Integer=10^4,x0::Array{Float64,1}=Float64[])
  const Mparams = M.params
  Mnoise! = M.noise!
  Mf! = M.f!
  function fn!(x::Array{Float64})
    Mf!(x,Mparams)
    Mnoise!(x)
    nothing
  end
  iterate(fn!,M,N,NI,x0)
end

function iterate(P::Peturbation,N::Int64=10^7,NI::Int64=10^4,Peps::Float64=P.defaulteps,x0::Array{Float64,1}=Float64[])
  const Pepsc = Peps
  const Mparams = P.M.params
  Mnoise! = P.M.noise!
  Mf! = P.M.f!
  function CX!(x::Array{Float64,1})
    x[:] = x + Pepsc*P.X(x)
    nothing
  end
  function fn!(x::Array{Float64,1})
    Mf!(x,Mparams)
    Mnoise!(x)
    CX!(x)
  end
  iterate(fn!,P.M,N,NI,x0)
end

# TODO: add params to noise fn

# Checking linear response

function findobservableslil(Peps::Float64,A::Array{Function,1},x::Array{Float64,1},P::Peturbation,N::Int64=10^7,NI::Int64=10^4)
  AN = length(A)
  xh = iterate(P,N,NI,Peps,x)
  x = xh[:,N]
  AV = Array(Float64,AN,2)
  for an = 1:AN
    AV[an,1] = mean(A[an](xh))
    AV[an,2] = observevar(A[an],xh)/length(xh)
  end
#  println("one obs done")
  return AV, x
end

function findobservablesbig(Peps::Float64,A::Array{Function,1},P::Peturbation,N::Int64=10^7,NI::Int64=10^4,NR::Int64=1)
  AN = length(A)
  AV = zeros(Float64,AN,2)
  x = Float64[]
  for r = 1:NR
    AV2,x = findobservableslil(Peps,A,x,P,N,NI)
    AV += AV2
  end
  AV[:,1] /= NR
  AV[:,2] /= NR^2
  return AV
end

function observepeturbation(P::Peturbation,A::Array{Function,1},epsv::Array{Float64,1},N::Int64=10^7,NI::Int64=10^4,NR::Int64=1)
  AN = length(A)
  eN = length(epsv)
  findobservablesbigf(Peps::Float64) = findobservablesbig(Peps,A,P,N,NI,NR)
  AX = pmap(findobservablesbigf,epsv)
#  println("iteration done")
  eA = Array(Float64,AN,eN)
  vA = Array(Float64,AN,eN)
  for i = 1:eN
    AXi = AX[i]
    eA[:,i] = AXi[:,1]
    vA[:,i] = AXi[:,2]
  end
  return eA,vA
end

function timedsample(P::Peturbation,PInitial::String,A::Array{Function,1},samplefn::Function,endtime::DateTime=tomorrowmorning(),N::Int64=10^7,NI::Int64=10^4,NR::Int64=1,NQ::Int64=3,samplefnargs=())
  AN = length(A)
  epsv = Array(Float64,0)
  eA = Array(Float64,AN,0)
  vA = Array(Float64,AN,0)
  NP = NQ * nworkers()
  filename = replace("RO-$(PInitial)$(endtime).h5",":","-")
  while now() < endtime
    epsvl = samplefn(NP,samplefnargs...)
    eAl, vAl = observepeturbation(P,A,epsvl,N,NI,NR)
    epsv = [epsv,epsvl]
    eA = [eA eAl]
    vA = [vA vAl]
    println("File saved at $(now())")
    save(filename,"epsv",epsv,"eA",eA,"vA",vA)
  end
  return epsv, eA, vA
end

function betasamplefn(NP::Int64,epsmax=0.0001,alpha=0.5,beta=1.)
  epsmax * rand(Distributions.Beta(alpha,beta),NP) # .*(1 - 2* (rand(Distributions.Bernoulli(),NP)))
end

function zerosamplefn(NP::Int64,epsmax=0.0001)
  zeros(Float64,NP)
end

function evensamplefn(NP::Int64,deps=0.0001)
  [0:NP-1]*deps
end

function checklinearresponse(epsv,eA,vA;epsmax=Inf,epsmin=-Inf,secondorder=false)
  (epsmax < Inf )|| (epsmin > -Inf) && ((epsv, eA, vA) = chopeps(epsv,eA,vA;epsmax=epsmax,epsmin=epsmin))
  AN = size(eA,1)
  eN = size(eA,2)
  errsize = sqrt(vA)
  neps = repmat(epsv',AN,1) ./ errsize
  neps2 = repmat(transpose(epsv.^2),AN,1) ./ errsize
  nA0c = 1 ./ errsize
  nAc = eA ./ errsize

  zeroval = Array(Float64,AN)
  lrtval = Array(Float64,AN)
  rss = Array(Float64,AN)
  pval = Array(Float64,AN)
  for an = 1:AN
    if secondorder
      lmX = [neps[an,:]' neps2[an,:]' nA0c[an,:]']
    else
      lmX = [neps[an,:]' nA0c[an,:]']
    end

    lmY = nAc[an,:]'
   # lmX |> println
    lmbh = inv(lmX' * lmX) * lmX' * lmY
    zeroval[an] = lmbh[1]
    lrtval[an] = lmbh[2]
    rss[an] = sum((lmY-lmX*lmbh).^2)
    pval[an] = Distributions.ccdf(Distributions.Chisq(eN-3),sum((lmY - lmX*lmbh).^2))
#    pval[an] = Distributions.ccdf(Distributions.Chisq(eN-2),sum((lmY - lmX*lmbh).^2))

#    println("Observable $an")
#    println("Coefficients: ", lmbh)
#    println("Residual sum of squares: ", sum((lmY-lmX*lmbh).^2))
#    println("p-value: ", Distributions.ccdf(Distributions.Chisq(eN-2),sum((lmY - lmX*lmbh).^2)))

#    lmX = [epsv ones(size(eA[an,:]'))]
#    lmY = eA[an,:]'
#    lmbh = inv(lmX' * lmX) * lmX' * lmY
#    println("Estimated variance given lr:", sqrt(sum((lmY-lmX*lmbh).^2)/(eN-2)))
  end
  return rss, pval, zeroval, lrtval
end

checklinearresponse(Dct::Dict,args...) = checklinearresponse(Dct["epsv"],Dct["eA"],Dct["vA"],args...)
