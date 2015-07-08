# Spectral determination of invariant measures

# Determines orbits of critical points of a map.
# Currently works only probably for 1 critical point
# and definitely only for 1 dimension.

function criticalorbit(M::IMap,Npts::Integer=300)
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
function logisticcriticalorbit(M::IMap,Npts::Integer=50)
  crit = BigFloat[0.5] #M.crit(M.params)
  alpha = BigFloat[M.params[1]]
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
                 dofirst::Bool=true, #include the first spike in the density
                 dorest::Bool=true) #include the rest of the spikes in the density
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

# Computes a spectral approximation of the acim
function spectralacim(M::IMap, # map whose acim we are finding
                      N::Integer=100; # number of points to take spectra at
                      verbose=false) # print output about accuracy etc
  crit = M.crit(M.params) #critical point(s), if any
  critexists = (length(crit) == 1) # is there a critical point?
  # (only going with one critical point for the moment)

  the_spectralpts = spectralpts(N,M.periodic,M.dom)
#  the_spectralpts[N] -= 2eps(1.) # otherwise get non-finite values on the boundary where there are spikes
#  the_spectralpts[1] += 2eps(1.)

  if critexists
    CO = criticalorbit(M)
    Sp = Spikes(CO,M.dom)

    eta_at_c = spikefn(crit,Sp)[1]
    Dr_at_c = spectralf(crit,[0:N-1],M.periodic,M.dom)

    zeta = inversetransfer(spikefn,M,(Sp,true,false))

    fixedhfn(x::Array{Float64,1}) =
      spikefn(x,Sp) - eta_at_c*zeta(x) # η - η(c)ζ
    h = transfer(fixedhfn,M,())(the_spectralpts) -
      spikefn(the_spectralpts,Sp,false,true) #L1(η - η(c)ζ) - η + η1
                                  # η1 is in here because later on we subtract the identity from
                                  # a matrix that looks like L1
  end

  if critexists
    fixedfn(x::Array{Float64,1},i::Int64) =
      spectralf(x,i,M.periodic,M.dom) - Dr_at_c[i+1]*zeta(x) # Ti - Ti(c)ζ (if we're using Chebyshev, e.g.)
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
    (M.dom[1] in Sp.CO.pts) && (LD[1,:] = 0; h[1] = 0)
    (M.dom[2] in Sp.CO.pts) && (LD[N,:] = 0; h[N] = 0)
  end

#  verbose && println("Doing linear algebra stuff")
  if critexists
    Lhat = [eta_at_c Dr_at_c; spectraltransf(N,M.periodic) * [h LD]]
  else
    Lhat = spectraltransf(N,M.periodic) * LD
  end

  Deltahat = [Lhat - I; zeros(M.Art.nfns,length(crit)) ac]
  weightm = speye(size(Deltahat,1)) # what the norm that you do the svd on looks like
  (U,S,V) = svd(weightm * Deltahat)
  verbose && println("Smallest singular values are ",[signif(S[i],4) for i = length(S) - [0:4]])
  r = V[:,end]


  # Creating the measure
  if critexists
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
    (M.Art.nfns > 0) && println("Artefact size: ",(Deltahat*r)[N+2])
  end

  return mu #, S # mu = measure, S = vector of singular values

end

