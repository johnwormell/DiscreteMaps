# Generic routines for measure and spectral functions

# Spectral helper functions - 1D ONLY
defdom(periodic) = periodic ? ([0 2pi]) : ([-1. 1.]) # default domain for Fourier/Chebyshev

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

function tridv(n::Integer,pow::Real=-1)
  tdv = zeros(n)
  tdv[2:2:end] = [1:floor(n/2)].^pow
  tdv
end

fourierint(n::Integer,dom::Array{Float64}=defdom(true)) = Tridiagonal(tridv(n-1,-1),zeros(n),-tridv(n-1,-1)) * domsize(dom)[1]/2pi
fourierdiff(n::Integer,dom::Array{Float64}=defdom(true)) = Tridiagonal(-tridv(n-1,1),zeros(n),tridv(n-1,1)) * 2pi/domsize(dom)[1]

# Fourier multiplication matrices
fillhf(m::Float64,f::Real) = fill(m,fld(f,2))

function fsmultk(k::Integer,N::Integer)
  # left,top
  I = [1,2k+1]
  J = [2k+1,1]
  V = [0.5,1]
  # top half
  append!(I,[2:2:N-2k-1,3:2:N-2k+1])
  append!(J,[2k+3:2:N,2k+2:2:N])
  append!(V,[fillhf(0.5,N-2k-1),fillhf(-0.5,N-2k)])

  # middle
  append!(I,[2:1:2k-1])
  append!(J,[2k-1:-1:2])
  append!(V,fill(0.5,2k-2))

  # bottom half
  append!(I,[2k+3:2:N,2k+2:2:N])
  append!(J,[2:2:N-2k-1,3:2:N-2k+1])
  append!(V,[fillhf(0.5,N-2k-1),fillhf(-0.5,N-2k)])

  return sparse(I,J,V,N,N)
end
function fcmultk(k::Integer,N::Integer)
  # left,top edges
  I = [1,2k]
  J = [2k,1]
  V = [0.5,1]
  # top half
  append!(I,[2:N-2k])
  append!(J,[2k+2:N])
  append!(V,fill(0.5,max(N-2k-1,0)))
  # middle
  append!(I,[2:2:2k-2,3:2:2k-1])
  append!(J,[2k-2:-2:2,2k-1:-2:3])
  append!(V,[fill(0.5,k-1),fill(-0.5,k-1)])
  # bottom half
  append!(I,[2k+2:N])
  append!(J,[2:N-2k])
  append!(V,fill(0.5,max(N-2k-1,0)))

  return sparse(I,J,V,N,N)
end
function fscmultk(k::Integer,N::Integer)
  (k==1)&& return speye(N)
  (rem(k,2) == 0) ? (return fcmultk(div(k,2),N)) : (return fsmultk(div(k-1,2),N))
end

function fscmult(coeffs::Array{Float64},N::Integer=length(coeffs))
  fscmultm = zeros(N,N)
  coeffsl = max(N,length(coeffs))
  for i in find(coeffs[1:coeffsl] .!= 0)
    fscmultm += coeffs[i] * fscmultk(i,N)
  end
  fscmultm
end

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
