# Spectral determination of invariant measures

# Determines orbits of critical points of a map.
# Currently works only probably for 1 critical point
# and definitely only for 1 dimension.

COdefaultNpts = 100

function criticalorbit(M::IMap,Npts::Integer=COdefaultNpts)
  crit = M.crit(M.params)
  critxx = M.critxx(M.params)

  Nc = length(crit)
  pts = Array(Float64,Npts,Nc)
  mag = Array(Float64,Npts,Nc)
  sgn = Array(Int8,Npts,Nc)

  for i = 1:Nc
    # the formulae here are described in Ruelle 2009, p. 1045
    pts[1,i] = M.f(crit[i],M.params)[1]
    mag[1,i] = sqrt(2/abs(critxx[i]))
    sgn[1,i] = sign(critxx[i])

    for p = 1:(Npts-1)
      pts[p+1,i:i] = M.f(pts[p,i],M.params)[1]
      dfp = (M.df(pts[p,i],M.params))[1]
      mag[p+1,i:i] = mag[p,i] ./ sqrt(abs(dfp))
      sgn[p+1,i:i] = sgn[p,i] * sign(dfp)
    end
  end
  # periodic error?
  return CriticalOrbit(crit,pts,mag,sgn)
end


# Using BigFloats to determine the critical orbit for logistic map -
# I think you need twice the precision or thereabouts for the positions
# of the orbits to have the same accuracy as the magnitudes - hence why this.
# Should probably let the domain of f, g etc just be Real to deal with this
# but it slows down iteration in the other parts of this module.
# Currently not being used for spectralacim.
function logisticcriticalorbit(M::IMap,Npts::Integer=COdefaultNpts)
  crit = BigFloat[M.crit(M.params)...] #
  alpha = BigFloat[M.params...]
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
# Computes a function containing all the spikes (ηi) from the map.
# The acim should be continuous after these spikes are removed.
function spikefn(x::Array{Float64,1},# points to evaluate at
                 Sp::Spikes, #spike measure
                 whichsp::Union(Integer,Array{Integer},Nothing)=nothing, # which spikes to do (nothing = all)
                 whichnsp::Union(Integer,Array{Integer},Nothing)=nothing, # if whichsp = nothing,
                 # which spikes not to do (nothing = all)
                 whichcp::Union(Integer,Array{Integer},Nothing)=nothing, # which cps to look at (nothing = all)
                 )
  # x = array of values to evaluate spike function at
  # Sp = container of spikes
  # dofirst = include first spikes (η1)
  # dorest = include the others spikes

  Nx = length(x)
  Nc = Sp.CO.Nc

  if whichsp == nothing
    spikeind = [1:Sp.CO.Npts]
    if whichnsp != nothing
      for i in whichnsp
        splice!(spikeind,i)
      end
    end
  else
    spikeind = [whichsp]
  end
  cpind = whichcp == nothing ? [1:Nc] : [whichcp]
  NSpts = length(spikeind)

  rawspikes = zeros(Float64,Nx)
  spikesum = zeros(Nx)

  for i in cpind
    for j = 1:NSpts
      rawspikes = Sp.mag0[i] * # size of initial bump
        abs(Sp.CO.mag[spikeind[j],i]) * # relative total size of spike vs bump
        (sign(x-Sp.CO.pts[spikeind[j],i]) .== Sp.CO.sgn[spikeind[j],i]) .*
        # is the x value on the right side of the spike?
        abs(x-Sp.CO.pts[spikeind[j],i]).^(-0.5) # size of spike at x value

      spikesum += indomain(x,Sp.dom) .*
        # is the x value in the domain of the map?
        rawspikes .*
        testfn(x,Sp.CO.pts[spikeind[j],i],Sp.widths[spikeind[j]])
        # test function centred around spikes to limit their support
    end
  end
  return spikesum
  #  dofull ? (return (full,first)) : (return first) # eek
end
spikefn(x::Array{Float64,1},Sp::Spikes;
        whichsp::Union(Integer,Nothing)=nothing,
        whichnsp::Union(Integer,Nothing)=nothing) = spikefn(x,Sp,whichsp,whichnsp)
spikefn(x::Float64,args...;kwargs...) = spikefn([x],args...;kwargs...)
measuredensity(x::F64U,Sp::Spikes) = spikefn(x,Sp)

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

# A right inverse of the transfer operator: used on η1 to calculate ζ (from the paper)
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

# Function to compute an invariant measure from a transfer operator
function findinvmeasure(Lhat::Array{Float64,2};verbose=false)
  Deltahat = copy(Lhat)
  for i = 1:minimum(size(Lhat))
    Deltahat[i,i] -= 1
  end

  weightm = speye(size(Deltahat,1)) # what the norm that you do the svd on looks like
  (U,S,V) = svd(weightm * Deltahat |> chopm)
  verbose && println("Smallest singular values are ",[signif(S[i],4) for i = length(S) - [0:4]])
  r = V[:,end]
  return (r,(U,S,V))
end


# Computes a spectral approximation of the acim
function spectralacim(M::IMap, # map whose acim we are finding
                      N::Integer=100; # number of points to take spectra at
                      verbose=false, # print output about accuracy etc
                      uselogisticcofn=false, # use logisticcriticalorbit function to calculate CO
                      returntransfmat=false, # if no critical points, additionally return the transfer operator matrix
                      sigma::Float64 = 0., # width of Gaussian noise
                      shortmatrix=(sigma==0), # for the spikes, return a collapsed matrix instead of a true one
                      usecrit=true, # use spikes/critical orbit stuff
                      lastspiketonoise = false, # turn the last spike into uniformly distributed noise
                      CONpts = COdefaultNpts) # number of points in critical orbit to use
#  M.Art = Artefacts()
  crit = usecrit ? M.crit(M.params) : [] #critical point(s), if any
  Nc = length(crit)
  critexists = (Nc > 0) # is there a critical point?
  noiseexists = (sigma > 0) # is there noise?
  if Nc > 1
    shortmatrix && error("Can't do collapsed critical points for more than one critical point")
    ~noiseexists && error("Can't do multiple critical points without noise")
  elseif Nc == 1
    ~noiseexists && ~shortmatrix && error("Can't do full matrix yet!")
  end

  the_spectralpts = spectralpts(N,M.periodic,M.dom)
  #  the_spectralpts[N] -= 2eps(1.) # otherwise get non-finite values on the boundary where there are spikes
  #  the_spectralpts[1] += 2eps(1.)

  if critexists
    noiseexists && (CONpts = 1)
    CO = uselogisticcofn ? logisticcriticalorbit(M,CONpts) : criticalorbit(M,CONpts)
    Sp = Spikes(CO,M.dom)

    zeta = Array(Function,Nc)
    for i = 1:Nc
      zeta[i] = inversetransfer(spikefn,M,(Sp,1,nothing,i))
    end
  end

  if critexists && ~noiseexists
    if shortmatrix
      eta_at_c = spikefn(crit,Sp)
    else
      eta_at_c = Array(Float64,Nc,Nc)
      for i = 1:Nc
        eta_at_c[:,i] = spikefn(crit,Sp,1,nothing,i)
      end
    end
#    eta_at_c = spikefn(crit,Sp)[1]

    fixedhfn = Array(Function,Nc)
    for i = 1:Nc
      tpts = lastspiketonoise ? (nothing, CONpts) : ()
      function fhfni(x::Array{Float64,1})
        zetav = Array(Float64,length(x),Nc)
        for k = 1:Nc
          for j = 1:length(x)
            zetav[j,k] = zeta[k]([x[j]])[1]
          end
        end

        spikefn(x,Sp,tpts...) - zetav * eta_at_c[:,i] # η - η(c)ζ
      end
      fixedhfn[i] = fhfni
    end

    if shortmatrix
      h = transfer(fixedhfn[1],M,())(the_spectralpts) -
        spikefn(the_spectralpts,Sp,nothing,1) #L1(η - η(c)ζ) - η + η1
      # η1 is in here because later on we subtract the identity from
      # a matrix that looks like L1
      if lastspiketonoise
        noisedist = ones(N) / sum(spectralvaluetotalint(N,M.periodic,M.dom))
        nconst = normalisedtestfnspiketotalintegral * Sp.mag0[1] *
         Sp.CO.mag[CONpts,1] .* sqrt(Sp.widths[CONpts,1])
        h += nconst * noisedist
      end
#     elseif noiseexists
#       h = Array(Float64,N,Nc)
#       for i = 1:Nc
#         h[:,i] = transfer(fixedhfn[i],M,())(the_spectralpts)
#       end
    end
  end

  if critexists
    Dr_at_c = spectralf(crit,[0:N-1],M.periodic,M.dom)

    function fixedfn(x::Array{Float64,1},i::Int64)
      cpremove = 0
      for j = 1:Nc
        cpremove += Dr_at_c[j,i+1] * zeta[j](x)
      end
      spectralf(x,i,M.periodic,M.dom) - cpremove # Ti - Ti(c)ζ (if we're using Chebyshev, e.g.)
    end
  else
    fixedfn(x::Array{Float64,1},i::Int64) =
      spectralf(x,i,M.periodic,M.dom)
  end

  LD = Array(Float64,N,N)
  ac = Array(Float64,M.Art.nfns,N)

  for i = 1:N
    LD[:,i] = transfer(fixedfn,M,(i-1))(the_spectralpts)

    # Getting the artefact functions and doing stuff with it
    if M.Art.nfns > 0
      ac[:,i] = M.Art.getcoeffs(
        transfer(fixedfn,M,(i-1))(M.Art.pointsin),
        i)
      LD[:,i] += M.Art.artefn(the_spectralpts) * ac[:,i]
    end
  end

  # floating point/nan error for cheby points on the boundary - for others we hope they don't exist
  if critexists
    (minabs(Sp.CO.pts - M.dom[1]) < 10*eps(M.dom[1])) && (LD[1,:] = 0; ~noiseexists && (h[1] = 0))
    (minabs(Sp.CO.pts - M.dom[2]) < 10*eps(M.dom[2])) && (LD[N,:] = 0; ~noiseexists && (h[N] = 0))
  end

  #  verbose && println("Doing linear algebra stuff")
  if critexists && ~noiseexists
    Lhat = [eta_at_c Dr_at_c; spectraltransf(N,M.periodic) * [h LD]; zeros(M.Art.nfns,length(crit)) ac]
  elseif critexists && noiseexists
    Lhat = [Dr_at_c; spectraltransf(N,M.periodic) * LD; ac]
  else
    Lhat = [spectraltransf(N,M.periodic) * LD; zeros(M.Art.nfns,length(crit)) ac]
  end

  chopm!(Lhat)

  if noiseexists
    spectralker = spectralgauskernel(N,M.periodic,M.dom,sigma) |> chopm
    critspikeker = Array(Float64,N,Nc)
    if critexists
      for i = 1:Nc
        critspikeker[:,i] = spikefn(the_spectralpts,Sp,1,nothing,i)
      end
      critspikeker = spectralker * (spectraltransf(N,M.periodic) * critspikeker)

#       for i = 1:Nc
#         dotfn = spectralvaluetotalint(N,M.periodic,M.dom) .* spikefn(the_spectralpts,Sp,1,nothing,i)
#         spikeint = normalisedtestfnspiketotalintegral*Sp.CO.mag[1,i] *Sp.mag0[i] * sqrt(Sp.widths[1,i])
#         for j = 1:N
#           smpol = spectralker[:,j]|>vec
#           smcoef = spectralapprox(crit[i],smpol,M.periodic,M.dom)[1]
#           smpol[1] -= smcoef
#           critspikeker[j,i] = dot(dotfn,
#                                   spectralapprox(the_spectralpts,smpol,M.periodic,M.dom)) +
#             smcoef * spikeint
#         end
#       end
    end
    if M.Art.nfns > 0
      artespikeker = spectralker * (spectraltransf(N,M.periodic) * M.Art.artefn(the_spectralpts))
    else
      artespikeker = Array(Float64,N,0)
    end
    totalker = [critspikeker spectralker|>full -artespikeker]
    Lhat = totalker * Lhat
  end

  (r,(U,S,V)) = findinvmeasure(Lhat,verbose=verbose)

  # Creating the measure
  if critexists && ~noiseexists
    mu = SumMeasure([SpectralMeasure(r[2:end],M.dom,M.periodic),r[1] * Sp])
  else
    mu = SpectralMeasure(r,M.dom,M.periodic)
  end

  # Normalising
  munorm = totalint(mu)[1]
  mu = mu/munorm

  # checking everything's ok
  if verbose
    println("r: ",round(r[1:10]'/munorm,4))
    println("Lr: ",round((Lhat*r)[1:10]'/munorm,4))
    ~noiseexists && (M.Art.nfns > 0) && println("Artefact size: ",(Lhat*r)[N+2])
  end

  returntransfmat && return (mu, Lhat)
  return mu #, S # mu = measure, S = vector of singular values

end

