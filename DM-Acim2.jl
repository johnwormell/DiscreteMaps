# Spectral helper functions - 1D ONLY
defdom(periodic) = periodic ? ([0 2pi]) : ([-1. 1.])

# Chebyshev functions and transforms - 1D ONLY. For descriptions of what they do see below

chebyinormcoords(x::F64U,dom::Array{Float64}) = (dom[1] + dom[2])/2 + x * (dom[2] - dom[1])/2  #  - 3*eps(dom[2]-dom[1])
chebynormcoords(x::F64U,dom::Array{Float64}) = (x - (dom[1] + dom[2])/2)/(dom[2] - dom[1]) * 2

chebyt(x::F64U,k::I64U,dom::Array{Float64}=defdom(false)) =cos(acos(chebynormcoords(x,dom)) * k')
chebypts(n::Integer,dom::Array{Float64}=defdom(false)) = chebyinormcoords(cos(linspace(-pi,0,n)),dom)
chebytransf(n::Integer) =
  [1,2*ones(n-1)] .* cos([0:(n-1)]*linspace(-pi,0,n)') .*
[1,2*ones(n-2),1]' / (2n-2)

chebyapprox(x::F64U,coeffs::F64U, dom::Array{Float64} = defdom(false)) =
  chebyt(x,[0:length(coeffs)-1],dom)*coeffs

chebytotalint(n::Integer,dom::Array{Float64,2} = defdom(false)) =
  [2., 0., -2 ./ ([2:n-1].^2 - 1)] .* (mod([0:n-1],2) .== 0) * domsize(dom)[1] / 2
chebyint(n::Integer, dom::Array{Float64,2} = defdom(false)) =
  Tridiagonal([2,1./[2:n-1]],zeros(n),[0,-1./[1:n-2]])/2  * domsize(dom)[1] / 2

# # Fourier functions and transforms - 1D ONLY - for explanations of what they do see below

fourierinormcoords(x::F64U,dom::Array{Float64}) = dom[1] + x * (dom[2] - dom[1])/2pi  #  - 3*eps(dom[2]-dom[1])
fouriernormcoords(x::F64U,dom::Array{Float64}) = (x - dom[1])/(dom[2] - dom[1]) * 2pi

function fourierscfn(k::Int64)
  (k == 0) && (return x-> 0*x+1)
  (rem(k,2) == 1) ? (return x -> cos((k+1)/2*x)) : (return x -> sin(k/2*x))
end

function fouriersc(x::F64U,k::I64U,dom::Array{Float64}=defdom(true))
  rmat = zeros(length(x),length(k))
  for i = 1:length(k)
    rmat[:,i] = fourierscfn(k[i])(fouriernormcoords(x,dom))
  end
  return rmat
end

fouriersc([0,pi,2pi],[0,1,2,3])
fourierpts(n::Integer,dom::Array{Float64}=defdom(true)) = fourierinormcoords([0:n-1]/n * 2pi,dom)
function fouriertransf(n::Integer)
  fp = fourierpts(n)
  ftm = fouriersc(fp,[0:(n-1)]) / n * 2
  ftm[:,1] /= 2
  return ftm'
end

fourierapprox(x::F64U,coeffs::F64U, dom::Array{Float64} = defdom(true)) =
  fouriersc(x,[0:length(coeffs)-1],dom)*coeffs

fouriertotalint(n::Integer,dom::Array{Float64}=defdom(true)) = [2pi, zeros(n-1)] * domsize(dom)[1] # check

function tridv(n::Integer)
  tdv = zeros(n)
  tdv[2:2:end] = 1./[1:floor(n/2)]
  tdv
end

fourierint(n::Integer,dom::Array{Float64}=defdom(true)) = Tridiagonal(tridv(n-1),zeros(n),-tridv(n-1)) * domsize(dom)[1]/2pi

# Generic spectral functions - 1D only

# Returns spectral functions of indices in k at values in x (e.g. Cheby T, sine/cosine)
spectralf(x::F64U,k::I64U,periodic=false,dom::Array{Float64}=defdom(periodic)) =
  periodic ? fouriersc(x,k, dom) : chebyt(x,k,dom)

# Returns sampling points of spectrum (even for Fourier, Cheby for Cheby)
spectralpts(n::Integer,periodic=false,dom::Array{Float64}=defdom(periodic)) =
  periodic ? fourierpts(n,dom) : chebypts(n,dom)
# Returns spectral transform matrix (i.e. turns value space into spectral space)
spectraltransf(n::Integer,periodic=false) =
  periodic ? fouriertransf(n) : chebytransf(n)

# Returns spectral approximation at values in x
spectralapprox(x::F64U,coeffs::F64U,periodic=false,dom::Array{Float64}=defdom(periodic)) =
  periodic ? (fourierapprox(x,coeffs,dom)) : (chebyapprox(x,coeffs,dom))

# Vector returns integral over domain when dotted with spectral coefficients
spectraltotalint(n::Integer,periodic=false,dom::Array{Float64}=defdom(periodic)) =
  periodic ? (fouriertotalint(n,dom)) : (chebytotalint(n,dom))

# Matrix turns coefficients of function into coefficients of antiderivative
spectralint(n::Integer,periodic = false,dom::Array{Float64}=defdom(periodic)) =
  periodic ? (fourierint(n,dom)) : (chebyint(n,dom))

# Measure object stuff

function measureintbd(intdom::Array{Float64,2},Sp::Spikes,A=nothing#, dofirst::Bool=true, dorest::Bool=true
                  )
  # x = array of values to evaluate spike function at
  # Sp = container of spikes
  # dofirst = include first spikes (η1)
  # dorest = include the others spikes
  Nc = Sp.CO.Nc
#   dorest && (rest = zeros(Nx))
#   dofirst && (first = zeros(Nx))

  theintegral = 0.
  theerror = 0.
  for i = 1:Nc
    lratio = restrictto(barrierdistoriented(Sp.CO.pts[:,i],Sp.CO.sgn[:,i],intdom[:,1]) ./ Sp.widths[:,i],0,1)
    rratio = restrictto(barrierdistoriented(Sp.CO.pts[:,i],Sp.CO.sgn[:,i],intdom[:,2]) ./ Sp.widths[:,i],0,1)
    partialspikes = (0 .< lratio .< 1) | (0 .< rratio .< 1)
    fullspikes = abs(lratio - rratio) .== 1

    quadgkspikes = find(partialspikes)
    # full spikes
    if A == nothing
      theintegral += normalisedtestfnspiketotalintegral * Sp.mag0[i] *
        sum(fullspikes .* Sp.CO.mag[:,i] .* sqrt(Sp.widths[:,i]))
    else
      append!(quadgkspikes,find(fullspikes))
    end

    # progressively add
    for j in quadgkspikes
      if A == nothing
        intfn = normalisedtestfnspike
      else
        intfn(x) = normalisedtestfnspike(x) * A((x - Sp.CO.pts[j,i])/Sp.widths[j,i] / Sp.CO.sgn[j,i])
      end
      (pint, perr) = quadgk(intfn,sort([lratio[j],rratio[j]])...)
      pcoeff = Sp.mag0[i] * Sp.CO.mag[j,i] .* sqrt(Sp.widths[j,i])
      theintegral += pint * pcoeff
      theerror += perr * pcoeff
    end
  end

  return theintegral, theerror
#   full  = zeros(Nx)
#   dofirst && (full += first)
#   dorest && (full += rest)
#   return full
  #  dofull ? (return (full,first)) : (return first) # eek
end

function measureintbd(intdom::Array{Float64,2},mu::SpectralMeasure,A=nothing)
  if A == nothing
    mint = diff(DiscreteMaps.spectralapprox(vec(intdom),DiscreteMaps.spectralint(mu.N,mu.periodic,mu.dom) * mu.coeffs,mu.periodic,mu.dom))[1]
    mu.periodic && (mint += mu.coeffs[1] * diff(vec(intdom))[1])
    return (mint,0.)
  else
    return quadgk(x -> A(x)*DiscreteMaps.spectralapprox(x,mu.coeffs,mu.periodic,mu.dom)[1],intdom[:]...)
  end
end

function measureintbd(intdom::Array{Float64,2},mu::SumMeasure,A=nothing)
  totalint = 0.
  totalerr = 0.
  for musub in mu.components
    mi = measureintbd(intdom,musub,A)
    totalint += mi[1]
    totalerr += mi[2]
  end
  return totalint, totalerr
end

function measureint(intdom::Array{Float64,1},mu::Measure,A=nothing)
  NB = length(intdom) - 1
  totalint = Array(Float64,NB)
  totalerr = Array(Float64,NB)
  for i = 1:NB
    (spi,spe) = DiscreteMaps.measureintbd([intdom[i] intdom[i+1]],mu,A)
    totalint[i] = spi
    totalerr[i] = spe
  end
  return (totalint, totalerr)
end

totalint(mu::Measure,A=nothing) = measureintbd(mu.dom,mu,A)
normalise(mu::Measure) = mu / totalint(mu)[1]
lyapunov(M::DifferentiableMap,mu::Measure) =
  exp(DiscreteMaps.totalint(mu,x->log(nextfloat(0.)+abs(M.df(x,M.params)[1])))[1]) # estimated Lyapunov exponent

# Measure density

# spikefunction measure is in DM-Acim.jl
measuredensity(x::F64U,mu::SpectralMeasure) = spectralapprox(x,mu.coeffs,mu.periodic,mu.dom)
function measuredensity(x::F64U,mu::SumMeasure)
  dens = zeros(size(x))
  for musub in mu.components
    dens += measuredensity(x)
  end
  return dens
end
