# Measure object stuff
zerosqrt(x::Array{Float64}) = sqrt(x .* (x .> 0))
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
  if Sp.widths == nothing
    if A == nothing
      theintegral += Sp.CO.mag .* Sp.mag0' .* Sp.CO.sgn .*
      2(zerosqrt(Sp.CO.sgn .* (intdom[2] - Sp.CO.pts)) - zerosqrt(Sp.CO.sgn .* (intdom[1] - Sp.CO.pts))) |> sum
    else
      intfn(x) = A(x)[1] * spikefn(x,Sp)[1]
      problempts = [intdom[1], sort(Sp.CO.pts[intdom[1] .< Sp.CO.pts .< intdom[2]]), intdom[2]]
      (pint, perr) = quadgk(intfn,problempts)
      theintegral += pint
      theerror += perr
    end
  else
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
          intfn(x) = normalisedtestfnspike(x) * A(x*Sp.widths[j,i]/Sp.CO.sgn[j,i] + Sp.CO.pts[j,i])
        end
        (pint, perr) = quadgk(intfn,sort([lratio[j],rratio[j]])...)
        pcoeff = Sp.mag0[i] * Sp.CO.mag[j,i] .* sqrt(Sp.widths[j,i])
        theintegral += pint * pcoeff
        theerror += perr * pcoeff
      end
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

function measureint(intdom::Array{Float64,1},mu::SpectralMeasure,A::Nothing)
  NB = length(intdom) - 1
  mint = diff(DiscreteMaps.spectralapprox(vec(intdom),DiscreteMaps.spectralint(mu.N,mu.periodic,mu.dom) * mu.coeffs,mu.periodic,mu.dom))
  mu.periodic && (mint += mu.coeffs[1] * diff(vec(intdom)))
  return (mint,zeros(NB))
end

function measureint(intdom::Array{Float64,1},mu::Measure,A=nothing)
  NB = length(intdom) - 1
  inttotal = Array(Float64,NB)
  errtotal = Array(Float64,NB)
  for i = 1:NB
    (spi,spe) = DiscreteMaps.measureintbd([intdom[i] intdom[i+1]],mu,A)
    inttotal[i] = spi
    errtotal[i] = spe
  end
  return (inttotal, errtotal)
end

totalint(mu::Measure,A::Function) = measureintbd(mu.dom,mu,A)
function totalint(mu::SumMeasure,A::Nothing=nothing)
  inttotal = 0.
  errtotal = 0.
  for musub in mu.components
    mi = totalint(musub,A)
    inttotal += mi[1]
    errtotal += mi[2]
  end
  return inttotal, errtotal
end
totalint(mu::SpectralMeasure,A::Nothing=nothing) = (dot(spectraltotalint(mu.N,mu.periodic,mu.dom),mu.coeffs),0.)
totalint(mu::Spikes,A::Nothing=nothing) = measureintbd(mu.dom,mu,A)

normalise(mu::Measure) = mu / totalint(mu)[1]
lyapunov(M::DifferentiableMap,mu::Measure) =
  exp(totalint(mu,x->log(nextfloat(0.)+abs(M.df(x,M.params)[1])))[1]) # estimated Lyapunov exponent

# Measure density

# spikefunction measure is in DM-Acim.jl
measuredensity(x::F64U,mu::SpectralMeasure) = spectralapprox(x,mu.coeffs,mu.periodic,mu.dom)
function measuredensity(x::F64U,mu::SumMeasure)
  dens = zeros(size(x))
  for musub in mu.components
    dens += measuredensity(x,musub)
  end
  return dens
end

# output is input for plotting a measure
function plotmeasure(mu::Measure;meshf=3000)
  boxsize = domsize(mu.dom)[1]/meshf
  hgd = linspace(mu.dom...,meshf+1)
  pgd = hgd[1:end-1] + boxsize/2
  mudens = measureint(hgd,mu)[1] / boxsize
  return pgd, mudens
end

# Fluctuation-Dissipation Theorem - calculating linear response of inv measure
function flucdis(mu::SpectralMeasure,L::Matrix{Float64},Xc::Array{Float64})#,Ac::Array{Float64})
  #  ~mu.periodic && error("Chebyshev measures not implemented for F-D")
  rhoc = mu.coeffs
  Ncoeffs = length(rhoc)
  divrX = spectraldiff(Ncoeffs,mu.periodic,mu.dom) * (spectralmult(Xc,Ncoeffs,mu.periodic)*rhoc)
  #  divrX2 = divrX[2:end]
  #  L2 = L[2:end,2:end]
  #   "max eigval: $(eigvals(L2) |> abs |> maximum)" |> println
  #   drho = divrX2
  #   for i = 1:4
  #     drho[:] = divrX2 + L2 * drho
  #   end
  drho = (I - L[2:end,2:end]) \ divrX[2:end]
  tint = spectraltotalint(Ncoeffs,mu.periodic,mu.dom)
  drho0 = -sum(tint[2:end] .* drho)/tint[1]
  return SpectralMeasure(-[drho0,drho],mu.dom,mu.periodic)
end
flucdis(mu::Measure,L::Matrix{Float64},X::Function) =
  flucdis(mu,L,spectraltransf(length(mu.coeffs),mu.periodic)
          *X(spectralpts(length(mu.coeffs),mu.periodic,mu.dom)))
function flucdis(M::IMap,XXc,N::Integer=100)
  (mu,L) = spectralacim(M,N,returntransfmat=true)
  flucdis(mu,L,XXc)
end
