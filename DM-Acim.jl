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
# Computes a function containing all the spikes (etai) from the map.
# The acim should be continuous after these spikes are removed.
function spikefn(x::Array{Float64,1},# points to evaluate at
                 Sp::Spikes, #spike measure
                 whichsp::Union(Integer,Array{Integer},Nothing)=nothing, # which spikes to do (nothing = all)
                 whichnsp::Union(Integer,Array{Integer},Nothing)=nothing, # if whichsp = nothing,
                 # which spikes not to do (nothing = all)
                 whichcp::Union(Integer,Array{Integer},Nothing)=nothing; # which cps to look at (nothing = all)
                 oneminustfn::Bool=false # use (1-phi) instead of phi (= bump fn) when multiplying the test function
                 )
  # x = array of values to evaluate spike function at
  # Sp = container of spikes
  # dofirst = include first spikes (eta1)
  # dorest = include the others spikes

  Nx = length(x)
  Nc = Sp.CO.Nc

  if whichsp == nothing
    spikeind = [1:Sp.CO.Npts]
    if whichnsp != nothing
      for i in sort([whichnsp],rev=true)
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

      tfnvals = vec(indomain(x',Sp.dom))
      # is the x value in the domain of the map?
      if Sp.widths != nothing
        tfnvals .*= testfn(x,Sp.CO.pts[spikeind[j],i],Sp.widths[spikeind[j]])
        # test function centred around spikes to limit their support
      end

      oneminustfn && (tfnvals = 1 - tfnvals)

      spikesum += tfnvals .* rawspikes
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

# A right inverse of the transfer operator: used on eta1 to calculate zeta (from the paper)
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
                      returntransfmat=false, # return the transfer operator matrix
                      returnsmallestsv=false, # return the smallest singular value
                      sigma::Float64 = 0., # width of Gaussian noise
                      shortmatrix=(sigma==0), # for the spikes, return a collapsed matrix instead of a true one
                      usecrit=true, # use spikes/critical orbit stuff
                      lastspiketonoise = false, # turn the last spike into uniformly distributed noise
                      CONpts = COdefaultNpts, # number of points in critical orbit to use
                      usetestfn = false # bound support of spikes with test functions
                      )
  crit = usecrit ? M.crit(M.params) : [] #critical point(s), if any
  Nc = length(crit)
  critexists = (Nc > 0) # is there a critical point?
  noiseexists = (sigma > 0) # is there noise?
  if Nc > 1
    shortmatrix && error("Can't do collapsed critical points for more than one critical point")
    ~noiseexists && error("Can't do multiple critical points without noise")
  elseif Nc == 1
    #   ~noiseexists && ~shortmatrix && error("Can't do full matrix yet!")
  end

  the_spectralpts = spectralpts(N,M.periodic,M.dom)
  #  the_spectralpts[N] -= 2eps(1.) # otherwise get non-finite values on the boundary where there are spikes
  #  the_spectralpts[1] += 2eps(1.)

  if critexists
    noiseexists && (CONpts = 1)
    CO = uselogisticcofn ? logisticcriticalorbit(M,CONpts) : criticalorbit(M,CONpts)
    Sp = usetestfn ? Spikes(CO,M.dom) : Spikes(CO,M.dom,fill(1.,Nc),nothing)#,min(CO.mag,0.4))

    zeta(x::Array{Float64,1},i::Integer) = inversetransfer(spikefn,M,(Sp,1,nothing,i))(x)
  end

  if critexists && ~noiseexists
    if shortmatrix
      tpts = lastspiketonoise ? (nothing, CONpts) : (nothing,nothing)
      eta_at_c = spikefn(crit,Sp,tpts...)
    else
      eta_at_c = Array(Float64,Nc,CONpts,Nc)
      for i = 1:Nc
        for j = 1:CONpts
          eta_at_c[:,j,i] = spikefn(crit,Sp,j,nothing,i)
        end
      end
    end
    #    eta_at_c = spikefn(crit,Sp)[1]

    if shortmatrix
      #        tpts = lastspiketonoise ? (nothing, CONpts) : (nothing,nothing)
      function fixedhfn(x::Array{Float64,1},i::Integer)
        i > Nc && error("Bounds error: i > Nc")
        zetav = Array(Float64,length(x),Nc)
        for k = 1:Nc
          for j = 1:length(x)
            zetav[j,k] = zeta([x[j]],k)[1]
          end
        end

        spikefn(x,Sp,tpts...,i) - zetav * eta_at_c[:,i] # eta - eta(c)zeta
      end
    else
      CONptslimit = lastspiketonoise ? CONpts-1 : CONpts
      function fixedhfn(x::Array{Float64,1},j::Integer,i::Integer)
        i > Nc && error("Bounds error: i > Nc")
        if (j == CONpts) && lastspiketonoise
          return zeros(length(x))
        elseif j > CONpts
          error("Bounds error: j > CONpts")
        end

        zetav = Array(Float64,length(x),Nc)
        for k = 1:Nc
          for l = 1:length(x)
            zetav[l,k] = zeta([x[l]],k)[1]
          end
        end

        spikefn(x,Sp,j,nothing,i) - zetav * eta_at_c[:,j,i] # eta - eta(c)zeta
      end
    end

    if shortmatrix
      h = transfer(fixedhfn,M,(1))(the_spectralpts) -
        spikefn(the_spectralpts,Sp,nothing,1) #L1(eta - eta(c)zeta) - eta + eta1
      # eta1 is in here because later on we subtract the identity from
      # a matrix that looks like L1
      if lastspiketonoise
        if Sp.widths == nothing
          nconst = 2sqrt(1-Sp.CO.pts[CONpts,1])
        else
          nconst = normalisedtestfnspiketotalintegral * sqrt(Sp.widths[CONpts,1])
        end

        h += fill(nconst * Sp.mag0[1] * Sp.CO.mag[CONpts,1] /
                    sum(spectralvaluetotalint(N,M.periodic,M.dom)),N)
      end
      #     elseif noiseexists
      #       h = Array(Float64,N,Nc)
      #       for i = 1:Nc
      #         h[:,i] = transfer(fixedhfn,M,(i))(the_spectralpts)
      #       end
    else
      h = Array(Float64,N,CONpts,Nc)
      for i = 1:Nc
        for j = 1:CONpts-1
          h[:,j,i] = #transfer(spikefn,M,(Sp,j,nothing,i))(the_spectralpts) -
          transfer(fixedhfn,M,(j,i))(the_spectralpts) -
            spikefn(the_spectralpts,Sp,j+1,nothing,i)
        end
        if lastspiketonoise
          if Sp.widths == nothing
            nconst = 2sqrt(1-Sp.CO.pts[CONpts,i])
          else
            nconst = normalisedtestfnspiketotalintegral * sqrt(Sp.widths[CONpts,i])
          end
          h[:,CONpts,i] += fill(nconst * Sp.mag0[i] * Sp.CO.mag[CONpts,i] /
                                  sum(spectralvaluetotalint(N,M.periodic,M.dom)),N)
        else
          h[:,CONpts,i] = transfer(fixedhfn,M,(CONpts,i))(the_spectralpts)
        end
      end
    end
  end

  if critexists
    Dr_at_c = spectralf(crit,[0:N-1],M.periodic,M.dom)

    function fixedfn(x::Array{Float64,1},i::Int64)
      cpremove = 0
      for j = 1:Nc
        cpremove += Dr_at_c[j,i+1] * zeta(x,j)
      end
      spectralf(x,i,M.periodic,M.dom) - cpremove # Ti - Ti(c)zeta (if we're using Chebyshev, e.g.)
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
  if critexists
    if noiseexists
      Lhat = [Dr_at_c; spectraltransf(N,M.periodic) * LD; ac]
    else
      if shortmatrix
        Lhat = [eta_at_c Dr_at_c; spectraltransf(N,M.periodic) * [h LD]; zeros(M.Art.nfns,length(crit)) ac]
      else
        topleftsize = Nc*CONpts
        topwidth = Nc*CONpts + N
        topm = zeros(topleftsize,topwidth)
        topm[CONpts*(0:Nc-1)+1,:] = [eta_at_c[:,:] Dr_at_c]
        for i = 1:Nc
          for j = 1:CONpts-1
            tlindex = CONpts*(i-1) + j
            topm[j+1,j] += 1
          end
        end

        Lhat = [topm; spectraltransf(N,M.periodic) * [h[:,:] LD]; zeros(M.Art.nfns,topleftsize) ac]
      end
    end
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
    if shortmatrix
      mu = SumMeasure([SpectralMeasure(r[2:end],M.dom,M.periodic),r[1] * Sp])
    else
      rsp = reshape(r[1:topleftsize],CONpts,Nc)
      sqrt(CONpts) * norm(std(rsp,2)) > 2S[end] * norm(r) && warn("Variance in spike coefficients is a bit large eh")
      Sp.mag0 .*= mean(rsp,1)|>vec
      mu = SumMeasure([SpectralMeasure(r[topleftsize+1:end],M.dom,M.periodic),Sp])
    end
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

  returntup = mu
  if returntransfmat
    if critexists & ~shortmatrix
      magconv = [vec(Sp.CO.mag),ones(N)]
      Lhat = (Lhat ./ magconv').*magconv
    end
    return returnsmallestsv ? (mu,Lhat,S[end]) : (mu,Lhat)
  else
    return returnsmallestsv ? (mu,S[end]) : mu
  end

end
