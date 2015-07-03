
# Determines orbits of critical points of a map.
# Currently works only probably for 1 critical point
# and definitely only for 1 dimension.

function criticalorbit(M::IMap,Npts::Integer=50)
  crit = BigFloat[0.5] #M.crit(M.params)
  alpha = BigFloat[3.8]
  critxx = M.critxx(M.params)

  Nc = length(crit)
  pts = Array(BigFloat,Npts,Nc)
  mag = Array(Float64,Npts,Nc)
  sgn = Array(Int8,Npts,Nc)

  for i = 1:Nc
    # the formulae here are described in Ruelle 2009, p. 1045
    pts[1,i] = (alpha .* crit .* (1 - crit))[1] #M.f(crit[i],M.params)[1]
    mag[1,i] = 1/sqrt(abs(critxx[i]/2))
    sgn[1,i] = sign(critxx[i])

    for p = 1:(Npts-1)
      pts[p+1,i:i] = alpha .* pts[p,i] .* (1 - pts[p,i]) #M.f(pts[p,i],M.params)
      dfp = (M.df(convert(Float64,pts[p,i]),M.params))
      mag[p+1,i:i] = mag[p,i] ./ sqrt(abs(dfp))
      sgn[p+1,i:i] = sgn[p,i] * sign(dfp)
    end
  end
  pts = convert(Array{Float64},pts)
  return CriticalOrbit(crit,pts,mag,sgn)
end

function logisticcospeeds(CO::CriticalOrbit,M::IMap)
  (Npts, Nc) = size(CO.pts)
  cospd = Array(Float64, Npts,Nc)
  crit = M.crit(M.params)
  cospd[1,1:1] = crit .* (1 - crit)
  for i = 1:(Npts-1)
    cospd[i+1,1:1] = 1 * (1 * CO.pts[i,:] .* (1-CO.pts[i,:])) +
      M.df(CO.pts[i,:],M.params) .* cospd[i,:]
  end
  return cospd
end

# Computes a function containing all the spikes from the map.
# The acim should be continuous after these spikes are removed.
function spikefn(x::Array{Float64,1},Sp::Spikes, dofirst::Bool=true, dorest::Bool=true)
  # x = array of values to evaluate spike function at
  # Sp = container of spikes
  # dofirst = include first spikes (η1)
  # dorest = include the others spikes
  Nx = length(x)
  Nc = Sp.CO.Nc
  xd = x'
  rawspikes = zeros(Float64,Sp.CO.Npts,Nx)
  dorest && (rest = zeros(Nx))
  dofirst && (first = zeros(Nx))
  for i = 1:Nc
    rawspikes[:,:] = (Sp.mag0[i] *
                        # size of initial bump
                        (sign(xd.-Sp.CO.pts[:,i]) .== Sp.CO.sgn[:,i]) .*
                      # is the x value on the right side of the spike?
                      abs(Sp.CO.mag[:,i]) # relative total size of spike vs bump
                      .* abs(xd.-Sp.CO.pts[:,i]).^(-0.5) # size of spike at x value
                      )
    dorest && (rest += (indomain(xd,Sp.dom) .*
                        # is the x value in the domain of the map?
                        sum(rawspikes[2:end,:] .*
                            testfn(xd,Sp.CO.pts[2:end,i],Sp.widths[2:end])
                            # test function centred around spikes to limit their support
                            ,1)) |> vec)

    dofirst && (first += indomain(xd,Sp.dom) .*
                rawspikes[1,:] .* testfn(xd,Sp.CO.pts[1,i],Sp.widths[1]) |> vec)

  end
  full  = zeros(Nx)
  dofirst && (full += first)
  dorest && (full += rest)
  return full
  #  dofull ? (return (full,first)) : (return first) # eek
end
spikefn(x::Array{Float64,1},Sp::Spikes; dorest::Bool=true, dofirst::Bool=true) = spikefn(x,Sp,dofirst,dorest)
spikefn(x::Float64,largs...) = spikefn([x],largs...)
measuredensity(x::F64U,Sp::Spikes) = spikefn(x,Sp)
#spikefn(x::Float64,Sp::Spikes,largs...) = spikefn([x],Sp,largs...)

# Transfer operator aka Frobenius-Perron operator
function transfer(r::Function,M::IMap,rargs=())
  function Lr(x::Array{Float64})
    Lrx = Array(Float64,size(x))
    for i = 1:length(x) #try using iterable stuff?????
      gxi = vec(M.g(x[i],M.params))
      Lrx[i] = sum(r(gxi,rargs...) ./
                   abs(M.df(gxi,M.params)))
    end

    return Lrx
  end
  Lr(x::Float64) = Lr([x])[1]
  return Lr
end

# A right inverse of the transfer
function inversetransfer(r::Function,M::IMap,rargs=())
  function Mr(x::Array{Float64})
    #    Mrx = Array(Float64,size(x))
    #     for i = 1:length(x)
    #       Mrx[i] = r(M.f(x[i],M.params),rargs...)[1] .* abs(M.df(x[i],M.params))[1] / branches
    #     end
    Mrx = r(M.f(x,M.params),rargs...).* abs(M.df(x,M.params)) .* M.inversetransferwt(x,M.params)
    return Mrx
  end
  return Mr
end

function spectralacim(M::IMap, # map whose acim we are finding
                      N::Integer=100; # number of points to take spectra at
                      verbose=true) # print output about accuracy etc
  crit = M.crit(M.params)
  crite = (length(crit) == 1) # only going with one critical point for the moment

  the_spectralpts = spectralpts(N,M.periodic,M.dom)
  the_spectralpts[N] -= 2eps(1.)
  the_spectralpts[1] += 2eps(1.)

  if crite
    CO = criticalorbit(M)
    Sp = Spikes(CO,M.dom)

    lphieta1 =  spikefn(crit,Sp,true,false)
    lphi = spectralf(crit,[0:N-1],M.periodic,M.dom)


    fixedhfn(x::Array{Float64,1}) =
      spikefn(x,Sp) -
      lphieta1[1]*inversetransfer(spikefn,M,(Sp,true,false))(x)
    h = transfer(fixedhfn,M,())(the_spectralpts) -
      spikefn(the_spectralpts,Sp,false,true)
  end

  if crite
    fixedfn(x::Array{Float64,1},i::Int64) =
      spectralf(x,i,M.periodic,M.dom) -
      lphi[i+1]*inversetransfer(spikefn,M,(Sp,true,false))(x)
  else
    fixedfn(x::Array{Float64,1},i::Int64) =
      spectralf(x,i,M.periodic,M.dom)
  end

  L = Array(Float64,N,N)
  ac = Array(Float64,M.Art.nfns,N)

  for i = 1:N
    L[:,i] = transfer(fixedfn,M,(i-1))(the_spectralpts)

    # Getting the artefact functions and doing stuff with it
    if M.Art.nfns > 0
      ac[:,i] = M.Art.getcoeffs(
        transfer(fixedfn,M,(i-1))(M.Art.pointsin),
        i)
      L[:,i] += M.Art.artefn(the_spectralpts) * ac[:,i]
    end
  end

  # floating point/nan error - to fix also for ~10 points away from boundary maybe
  if crite
    (M.dom[1] in Sp.CO.pts) && (L[1,:] = 0; h[1] = 0)
    (M.dom[2] in Sp.CO.pts) && (L[N,:] = 0; h[N] = 0)
  end

  verbose && println("Doing linear algebra stuff")
  if crite
    LL = [lphieta1 lphi; spectraltransf(N,M.periodic) * [h L]]
  else
    LL = spectraltransf(N,M.periodic) * L
  end

  LIF = [LL - I; zeros(M.Art.nfns,length(crit)) ac]
  weightv = ones(size(LIF)) #max([ones(length(crit)),1:N,ones(M.Art.nfns)],15) #max([-1:N-1],10-[-1:N-1])
  (U,S,V) = svd(weightv .* LIF)
  r = V[:,end]
  verbose && println("Smallest singular values are ",[signif(S[i],4) for i = length(S) - [0:4]])

  # Creating the measure
  if crite
    mu = SumMeasure([SpectralMeasure(r[2:end],M.dom,M.periodic),r[1] * Sp])
  else
    mu = SpectralMeasure(r,M.dom,M.periodic)
  end

  # Normalising
  munorm = totalint(mu)[1]
  mu = mu/munorm

  #  crite && (Sp = Spikes(Sp,normfactor=1. /rnorm) ) <- not this:
  #                                  # for the moment, you always put r[1] in front of spikes

  if verbose
    println("r: ",round(r[1:10]'/munorm,4))
    println("Lr: ",round((LL*r)[1:10]'/munorm,4))
    (M.Art.nfns > 0) && println("Artefact size: ",(LIF*r)[N+2])
  end

  return mu, S

end

