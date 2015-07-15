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

function chebydiff(n::Integer, dom::Array{Float64,2} = defdom(false))
  diffm = zeros(n,n)
  for i = 1:n
    if rem(i,2) == 0
      diffm[1,i] = i-1
      diffm[3:2:i-1,i] = 2i-2
    else
      diffm[2:2:i-1,i] = 2i-2
    end
  end
  diffm * 2/domsize(dom)[1]
end

function chebymult(coeffs::Array{Float64}, n::Integer=length(coeffs))
  coeffs2 = [coeffs[2:end]/2,zeros(max(0,n-length(coeffs)))] # length n-1
  coeffstpl = [coeffs2|>flipud,coeffs[1],coeffs2] # length 2n-1

  multm = zeros(n,n)
  for i = 1:n
    multm[1,i] += coeffstpl[n-i+1]
    for j = 2:n-i+1
      multm[j,i] += coeffstpl[n-i+j] + coeffs2[i+j-2]
    end
    for j = n-i+2:n
        multm[j,i] += coeffstpl[n-i+j]
    end
  end
  multm
end

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
function fouriersmultk(k::Integer,n::Integer)
  # left,top
  I = [1,2k+1]
  J = [2k+1,1]
  V = [0.5,1]
  # top half
  append!(I,[2:2:n-2k-1,3:2:n-2k+1])
  append!(J,[2k+3:2:n,2k+2:2:n])
  append!(V,[fillhf(0.5,n-2k-1),fillhf(-0.5,n-2k)])

  # middle
  append!(I,[2:1:2k-1])
  append!(J,[2k-1:-1:2])
  append!(V,fill(0.5,2k-2))

  # bottom half
  append!(I,[2k+3:2:n,2k+2:2:n])
  append!(J,[2:2:n-2k-1,3:2:n-2k+1])
  append!(V,[fillhf(0.5,n-2k-1),fillhf(-0.5,n-2k)])

  return sparse(I,J,V,n,n)
end
function fouriercmultk(k::Integer,n::Integer)
  # left,top edges
  I = [1,2k]
  J = [2k,1]
  V = [0.5,1]
  # top half
  append!(I,[2:n-2k])
  append!(J,[2k+2:n])
  append!(V,fill(0.5,max(n-2k-1,0)))
  # middle
  append!(I,[2:2:2k-2,3:2:2k-1])
  append!(J,[2k-2:-2:2,2k-1:-2:3])
  append!(V,[fill(0.5,k-1),fill(-0.5,k-1)])
  # bottom half
  append!(I,[2k+2:n])
  append!(J,[2:n-2k])
  append!(V,fill(0.5,max(n-2k-1,0)))

  return sparse(I,J,V,n,n)
end
function fouriermultk(k::Integer,n::Integer)
  (k==1)&& return speye(n)
  (rem(k,2) == 0) ? (return fouriercmultk(div(k,2),n)) : (return fouriersmultk(div(k-1,2),n))
end

function fouriermult(coeffs::F64U,n::Integer=length(coeffs))
  multm = zeros(n,n)
  coeffsl = min(n,length(coeffs))
  for i in find(coeffs[1:coeffsl] .!= 0)
    multm += coeffs[i] * fouriermultk(i,n)
  end
  multm
end

function fourierconv(coeffs::Array{Float64,1},dom::Array{Float64,2}, n::Integer=length(coeffs))
  coeffsl = min(n,length(coeffs))
  coeffsf = coeffs[1:coeffsl]
  coeffsl < n && append!(coeffsf,zeros(n - coeffsl))

  d = Array(Float64,n)
  d[1] = coeffsf[1]
  for i = 2:n
    d[i] = coeffsf[2fld(i,2)]/2
    end
  dl = Array(Float64,n-1)
  for i = 1:n-1
    dl[i] = rem(i,2) == 0 ? coeffsf[2div(i,2)+1]/2 : 0
  end

  Tridiagonal(dl,d,-dl)*domsize(dom)[1]
end

# Generic spectral functions - 1D only

# Returns spectral functions of indices in k at values in x (e.g. Cheby T, sine/cosine)
spectralf(x::F64U,k::I64U,periodic::Bool=false,dom::Array{Float64}=defdom(periodic)) =
  periodic ? fouriersc(x,k, dom) : chebyt(x,k,dom)

# Returns sampling points of spectrum (even for Fourier, Cheby for Cheby)
spectralpts(n::Integer,periodic::Bool=false,dom::Array{Float64}=defdom(periodic)) =
  periodic ? fourierpts(n,dom) : chebypts(n,dom)
# Returns spectral transform matrix (i.e. turns value space into spectral space)
spectraltransf(n::Integer,periodic::Bool=false) =
  periodic ? fouriertransf(n) : chebytransf(n)

# Returns spectral approximation at values in x
spectralapprox(x::F64U,coeffs::F64U,periodic::Bool=false,dom::Array{Float64}=defdom(periodic)) =
  periodic ? (fourierapprox(x,coeffs,dom)) : (chebyapprox(x,coeffs,dom))

# Vector returns integral over domain when dotted with spectral coefficients
spectraltotalint(n::Integer,periodic::Bool=false,dom::Array{Float64}=defdom(periodic)) =
  periodic ? (fouriertotalint(n,dom)) : (chebytotalint(n,dom))

# Matrix turns coefficients of function into coefficients of antiderivative
spectralint(n::Integer,periodic::Bool = false,dom::Array{Float64}=defdom(periodic)) =
  periodic ? (fourierint(n,dom)) : (chebyint(n,dom))

# Matrix turns coefficients of function into coefficients of derivative
spectraldiff(n::Integer,periodic::Bool = false,dom::Array{Float64}=defdom(periodic)) =
  periodic ? (fourierdiff(n,dom)) : (chebydiff(n,dom))

# Matrix turns coefficients of function into coefficients of function
# multiplied by a function approximated by coeffs
spectralmult(coeffs::Array{Float64},n=length(coeffs),periodic::Bool = false) =
  periodic ? (fouriermult(coeffs,n)) : (chebymult(coeffs,n))

# Fourier coefficients of periodic Gaussian distribution
defaultrelsigmasize = 0.001
function gaussianfouriercoefs(n::Integer,dom::Array{Float64,2}=DM.defdom(true),
                              sigma::Float64=defaultrelsigmasize*DM.domsize(dom)[1],mu::Float64=0.)
  ds = DM.domsize(dom)[1]
  omega0 = mu*2pi/ds
  gfv = Array(Float64,n)
  gfv[1] = 1. / ds
  for i = 2:n
    i2 = fld(i,2)
    gfv[i] = (rem(i,2)==0 ? cos(-i2*omega0) : sin(-i2*omega0)) *
      exp(-(2pi/ds * sigma * i2)^2 / 2) * (2 / ds)
  end
  gfv
end
function chebygauskernel(n::Integer, sratio::Float64=0.001)
#   omega0 = pi/2
# #  -(2pi*sratio*[0:n-1]).^2 / 2 |> exp |> println
#   kerm = diagm(cos(-[0:n-1]*omega0).*exp(-(2pi*sratio*[0:n-1]).^2 / 2))
#   kerm
  cm = chebytransf(n)
  cp = chebypts(n,[-1. 1.])
  km = Array(Float64,n,n)
  for i = 1:n
    for j = 1:n
      km[j,i] = (1 + (j!=1))/pi * sqrt(1 - cp[j]^2) * (exp(-((cp[j]-cp[i]) / 2sratio).^2 / 2) +
                   exp(-((cp[j]+cp[i]-2) / 2sratio).^2 / 2) +
                   exp(-((cp[j]+cp[i]+2) / 2sratio).^2 / 2))/sqrt(2pi)/2sratio
    end
  end
  return km #cm * km * cm'
end
# Returns kernel matrix of a Fourier Gaussian distribution
fouriergauskernel(n::Integer,sratio::Float64=0.001) =
    DM.fourierconv(gaussianfouriercoefs(n,[0. 2pi],2pi*sratio,0.),[0. 2pi],n)
fouriergauskernel(n::Integer,dom::Array{Float64,2} = DM.defdom(true),
                   sigma::Float64=defaultrelsigmasize*DM.domsize(dom)[1]) =
  DM.fourierconv(gaussianfouriercoefs(n,dom,sigma,0.),[0. 2pi],n)
