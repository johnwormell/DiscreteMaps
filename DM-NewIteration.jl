function oiterationpreloop!(x::Array{Float64},fn!::Function,NI::Integer)
  for n = 1:NI
    fn!(x)
  end
  nothing
end

function oiterationhistoryloop!(x_history::Array{Float64},x::Array{Float64},fn!::Function,NH::Integer)
  for n = 1:NH
    fn!(x)
    x_history[:,n] = x
  end
  nothing
end

# function niterationacvloop!(x::Array{Float64},fn!::Function,A::Array{Function,1},sA::Array{Float64,1},sAA::Array{Float64,2},N::Integer)
# # sA = (1x)AN, sAA = NL x AN
#   size(sAA,1) = NL

#   for n = 1:N
#     fn!(x)
#     for an = 1:AN
#       sA[an] += A[an](x)
#     end
#     x, sA
#   end


function oiterationmainloop!(x::Array{Float64},fn!::Function,A::Array{Function,1},sA::Array{Float64,1},N::Integer)
  AN = length(A)
  for n = 1:N
    fn!(x)
    for an = 1:AN
      sA[an] += A[an](x)[1]
    end
  end
  nothing
end

function obsiterate(fn!::Function,It::IterationSchema;returnhistory=false)
  NH = min(It.NH,It.N)
  x_history = Array(Float64,It.P.M.dim,NH)
  x = It.P.M.init()
  oiterationpreloop!(x,fn!,It.NI)

  oiterationhistoryloop!(x_history,x,fn!,NH)
  returnhistory && return x_history

  sA = Array(Float64,It.AN)
  vA = Array(Float64,It.AN)
  for an = 1:It.AN
    sA[an] = sum(It.A[an](x_history))
    vA[an] = observevar(It.A[an],x_history)/It.N
  end
  NH < It.N && oiterationmainloop!(x,fn!,It.A,sA,It.N-NH)
  eA = sA/It.N
  return [eA vA]
end

function obsiterate(It::IterationSchema, Peps::Float64=0.;returnhistory=false)
  const Pepsc = Peps
  const Mparams = It.P.M.params
  Mnoise! = It.P.M.noise!
  Mf! = It.P.M.f!
  function CX!(x::Array{Float64,1})
    x[:] = x + Pepsc*It.P.X(x)
    nothing
  end
  function fn!(x::Array{Float64,1})
    Mf!(x,Mparams)
    CX!(x) # swapped this with Mnoise!
    Mnoise!(x)
  end
  obsiterate(fn!,It;returnhistory=returnhistory)
end

function timedsample(It::IterationSchema; endtime::DateTime=tomorrowmorning(),
                     NQ::Int64=1, NP::Int64=NQ*nworkers(), NCycles::Real=Inf,
                     startstring::String="rs")
  epsv = Array(Float64,0)
  eAv = Array(Float64,It.AN,0)
  vAv = Array(Float64,It.AN,0)
  filename = replace("$(startstring)-$(It.PInitial)-$(endtime)--$(hour(now()))-$(minute(now())).h5",":","-")
  obsiteratef(Peps::Float64) = obsiterate(It,Peps)
  cyclecount = 0
  while (now() < endtime) & (cyclecount < NCycles)
    epsl = It.samplefn(NP,It.samplefnargs...)
    pmapcontainer = pmap(obsiteratef,epsl)
    eAv = [eAv Array(Float64,It.AN,NP)]
    vAv = [vAv Array(Float64,It.AN,NP)]
    epsv = [epsv, epsl]
    for i = 1:NP
      dA = pmapcontainer[i]
      eAv[:,cyclecount*NP+i] = dA[:,1]
      vAv[:,cyclecount*NP+i] = dA[:,2]
    end

    save(filename,"epsv",epsv,"eA",eAv,"vA",vAv,"N",It.N,"NH",It.NH)
    println("File saved at $(now())")
    cyclecount += 1
  end
  return cyclecount, epsv, eAv, vAv
end
timedsample(PInitial::String; endtime::DateTime=tomorrowmorning(),
            NQ::Int64=1, NP::Int64=NQ*nworkers(), NCycles::Real=Inf,
            Itargs=(),Itkwargs=kw(),
            startstring::String="rs") =
  timedsample(itdict[PInitial](Itargs...;Itkwargs...);endtime=endtime,NQ=NQ,NP=NP,NCycles=NCycles,
              startstring=startstring)
