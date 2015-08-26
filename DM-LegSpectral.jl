# Generic routines for measure and spectral functions
# Packages
try
  Pkg.installed("FastGaussQuadrature")
catch
  error("Maybe you need to install FastGaussQuadrature, available at https://github.com/ajt60gaibb/FastGaussQuadrature.jl.")
end
using FastGaussQuadrature

# Spectral helper functions - 1D ONLY
defdom(periodic) = periodic ? ([0 2pi]) : ([-1. 1.]) # default domain for Fourier/Legendre

# Legendre functions and transforms - 1D ONLY. For descriptions of what they do see below

leginormcoords(x::F64U,dom::Array{Float64}) = (dom[1] + dom[2])/2 + x * (dom[2] - dom[1])/2  #  - 3*eps(dom[2]-dom[1])
legnormcoords(x::F64U,dom::Array{Float64}) = (x - (dom[1] + dom[2])/2)/(dom[2] - dom[1]) * 2

function legp(x::F64U,k::I64U,dom::Array{Float64}=defdom(false))

  xn = legnormcoords(x,dom)
  ks = issorted(k) ? k : sort(k)

  lpmat = Array(Float64,length(xn),length(k))
  maxk = ks[end]
  maxkind = length(k)
  kind = 1

  p2 = ones(length(xn))
  while kind <= maxkind && ks[kind] == 0
    lpmat[:,kind] = p2
    kind += 1
  end

  p1 = [x]
  while kind <= maxkind && ks[kind] == 1
    lpmat[:,kind] = p1
    kind += 1
  end

  p0 = Array(Float64,length(xn))
  for j = 2:maxk
    for i = 1:length(xn)
      p0[i] = ((2j-1)*xn[i]*p1[i] - (j-1)*p2[i])/j
    end

    while kind <= maxkind && ks[kind] == j
      lpmat[:,kind] = p0
      kind += 1
    end

    p2[:] = p1[:]
    p1[:] = p0[:]
  end

  issorted(k) || (lpmat = lpmat[:,sortperm(k)])

  lpmat
end

legpts(n::Integer,dom::Array{Float64}=defdom(false)) = leginormcoords(FastGaussQuadrature.gausslegendre(n)[1],dom)
function legtransf(n::Integer)
  (x, w) = FastGaussQuadrature.gausslegendre(n)
  legp(x,[0:n-1])'.*[0.5:1:n-0.5].*w'
end

legapprox(x::F64U,coeffs::F64U,dom::Array{Float64} = defdom(false)) =
  legp(x,[0:length(coeffs)-1],dom) * coeffs

legtotalint(n::Integer,dom::Array{Float64,2} = defdom(false)) =
  [1, zeros(n-1)] * domsize(dom)[1]
legvaluetotalint(n::Integer,dom::Array{Float64,2} = defdom(false)) =
  FastGaussQuadrature.gausslegendre(n)[2] * domsize(dom)[1]/2

leginnerprodm(n::Integer=length(coeffs),dom::Array{Float64,2}=defdom(false)) =
  Diagonal(domsize(dom)[1]./[1:2:2n-1])

legint(n::Integer, dom::Array{Float64,2} = defdom(false)) =
  Tridiagonal(1./[1:2:2n-3],zeros(n),-1./[1:2:2n-3]) * domsize(dom)[1] / 2

function legdiff(n::Integer, dom::Array{Float64,2} = defdom(false))
  diffm = zeros(n,n)
  for i = 1:n
    if rem(i,2) == 0
      diffm[1:2:i-1,i] = [1:4:2i-3]
    else
      diffm[2:2:i-1,i] = [3:4:2i-3]
    end
  end
  diffm * 2/domsize(dom)[1]
end

function legmult(coeffs::Array{Float64}, n::Integer=length(coeffs))
#   coeffs2 = [coeffs[2:end]/2,zeros(max(0,n-length(coeffs)))] # length n-1
#   coeffstpl = [coeffs2|>flipud,coeffs[1],coeffs2] # length 2n-1

  mconstarray2 = eye(n)
  multm = coeffs[1] * mconstarray2
  mconstarray1 = Tridiagonal(0.5+1./[2:4:4n-6],zeros(n),(0.5-1./[6:4:4n-2])) |> full
  coeffs[2] != 0 && (multm += coeffs[2] * mconstarray1)
  for j = 2:maximum([find(coeffs.!=0),2])-1
    mconstarray = -mconstarray2*(j-1)/j

    for i = 1:n-1
      mconstarray[:,i+1] += (2j-1)/j/(2i+1) * i*mconstarray1[:,i]
    end
    for i = 0:n-2
      mconstarray[:,i+1] += (2j-1)/j/(2i+1) * (i+1)*mconstarray1[:,i+2]
    end

    coeffs[j+1] != 0 && (multm += coeffs[j+1] * mconstarray)
    mconstarray2 = mconstarray1
    mconstarray1 = mconstarray
  end
  multm
end

function legconvgetptmatrix(coeffs::Array{Float64,1},n::Integer=length(coeffs))
  bleft = Array(Float64,n,n)
  bleft[1,1] = coeffs[1] - coeffs[2]/3
  for i = 1:n-2
    bleft[i+1,1] = coeffs[i]/(2i-1) - coeffs[i+2]/(2i+3)
  end
  bleft[n,1] = coeffs[n-1]/(2n-3)

  bleft[1,2] = - bleft[2,1]/3
  for i = 1:n-2
    bleft[i+1,2] = bleft[i,1]/(2i-1) - bleft[i+1,1] - bleft[i+2,1]/(2i+3)
  end
  bleft[n,2] = bleft[n-1,1]/(2n-3) - bleft[n,1]

  for j = 2:n-1
    for i = 0:j-1
      bleft[i+1,j+1] = (-1)^(i+j) * (2i+1)/(2j+1) * bleft[j+1,i+1]
    end

    for i = j:n-2
      bleft[i+1,j+1] = (2j-1)/(2i-1) * bleft[i,j] -(2j-1)/(2i+3) * bleft[i+2,j] + bleft[i+1,j-1]
    end
    bleft[n,j+1] = (2j-1)/(2n-3) * bleft[n-1,j] + bleft[n,j]
  end
  bleft
end

function legconv(coeffs::Array{Float64,1},dom::Array{Float64,2}=defdom(false), n::Integer=length(coeffs))
  bleft = legconvgetptmatrix(coeffs,n)
  bright = legconvgetptmatrix(coeffs.*(-1).^[0:n-1],n).*(-1).^[0:n-1]'
  lpts = legpts(n)
  ptbrk = div(n,2)

  legtransf(n)*
    [legp(lpts[1:ptbrk]+1,[0:2:n-1])*bleft[1:2:n,:],
     legp(lpts[ptbrk+1:end]-1,[0:2:n-1])*bright[1:2:n,:]]*domsize(dom)[1]
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

fouriertotalint(n::Integer,dom::Array{Float64}=defdom(true)) = [1, zeros(n-1)] * domsize(dom)[1]
fouriervaluetotalint(n::Integer,dom::Array{Float64}=defdom(true)) = fill(domsize(dom)[1]/n,n)

function fourierinnerprodm(n::Integer=length(coeffs),dom::Array{Float64,2}=defdom(true))
  fipd = Diagonal(fill(domsize(dom)[1]/2,n))
  fipd.diag[1] *= 2
  fipd
end

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

function fourierconv(coeffs::Array{Float64,1},dom::Array{Float64,2}=defdom(true), n::Integer=length(coeffs))
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
  periodic ? fouriersc(x,k, dom) : legp(x,k,dom)

# Returns sampling points of spectrum (even for Fourier, Gauss-Legendre nodes for Legendre)
spectralpts(n::Integer,periodic::Bool=false,dom::Array{Float64}=defdom(periodic)) =
  periodic ? fourierpts(n,dom) : legpts(n,dom)
# Returns spectral transform matrix (i.e. turns value space into spectral space)
spectraltransf(n::Integer,periodic::Bool=false) =
  periodic ? fouriertransf(n) : legtransf(n)

# Returns spectral approximation at values in x
spectralapprox(x::F64U,coeffs::F64U,periodic::Bool=false,dom::Array{Float64}=defdom(periodic)) =
  periodic ? (fourierapprox(x,coeffs,dom)) : (legapprox(x,coeffs,dom))

# Vector returns integral over domain when dotted with spectral coefficients
spectraltotalint(n::Integer,periodic::Bool=false,dom::Array{Float64}=defdom(periodic)) =
  periodic ? (fouriertotalint(n,dom)) : (legtotalint(n,dom))

# Vector returns integral over domain when dotted with values at spectral points
spectralvaluetotalint(n::Integer,periodic::Bool=false,dom::Array{Float64}=defdom(periodic)) =
  periodic ? (fouriervaluetotalint(n,dom)) : (legvaluetotalint(n,dom))

# Vector returns L2 inner product matrix for spectral coefficients
spectralinnerprodm(n::Integer,periodic::Bool=false,dom::Array{Float64}=defdom(periodic)) =
  periodic ? (fourierinnerprodm(n,dom)) : (leginnerprodm(n,dom))

# Matrix turns coefficients of function into coefficients of antiderivative
spectralint(n::Integer,periodic::Bool = false,dom::Array{Float64}=defdom(periodic)) =
  periodic ? (fourierint(n,dom)) : (legint(n,dom))

# Matrix turns coefficients of function into coefficients of derivative
spectraldiff(n::Integer,periodic::Bool = false,dom::Array{Float64}=defdom(periodic)) =
  periodic ? (fourierdiff(n,dom)) : (legdiff(n,dom))

# Matrix turns coefficients of function into coefficients of function
# multiplied by a function approximated by coeffs
spectralmult(coeffs::Array{Float64},n=length(coeffs),periodic::Bool = false) =
  periodic ? (fouriermult(coeffs,n)) : (legmult(coeffs,n))

# Matrix turns coefficients of function into coefficients of function
# convolved with a function approximated by coeffs
spectralconv(coeffs::Array{Float64},periodic::Bool = false,dom::Array{Float64}=defdom(periodic),n=length(coeffs)) =
  periodic ? (fourierconv(coeffs,dom,n)) : (legconv(coeffs,dom,n))


# Fourier coefficients of periodic Gaussian distribution
defaultrelsigmasize = 0.001
function fouriergauscoefs(n::Integer,dom::Array{Float64,2}=defdom(true),
                              sigma::Float64=defaultrelsigmasize*domsize(dom)[1],mu::Float64=0.)
  ds = domsize(dom)[1]
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

legsigmaratbuf = 1.5
legsigmaratwarn = 6.

# Legendre coefficients of scaled chopped Gaussian distribution
function leggauscoefs(n::Integer,dom::Array{Float64,2}=defdom(false),
                              sigma::Float64=defaultrelsigmasize*domsize(dom)[1],mu::Float64=mean(dom))
  legnormcoords(mu,dom) == 0. || error("leggauscoeffs does not support non-centered mu values")
  ds = domsize(dom)[1]
  sigmarat = sigma * 2 / domsize(dom)[1]
  sigmarat < legsigmaratwarn/n && warn("The noise size may be too small for the number of spectral coefficients")
  if sigmarat >= sqrt(-1/2log(eps(1.)))/legsigmaratbuf
    # for a wide peak we can do with a standard Legendre transform
    gfv = legtransf(n)*exp(-(legpts(n,dom)-mu).^2/2sigma^2)/(sigma*sqrt(2pi))
  else
    # for a narrow peak we only look at the effective support of the gaussian to get a better transform
    normwdth = sigmarat * sqrt(-log(eps(1.)))
    pkn = int(n / normwdth)
    (pkpts,pkw) = FastGaussQuadrature.gausslegendre(pkn)
    pknrange = int((pkn-n)/2)+(1:n)
    pkpts = pkpts[pknrange]
    pkw = pkw[pknrange]
    gfv = exp(-pkpts.^2/2sigmarat^2)/(sigmarat*sqrt(2pi))
    gfv = (legp(pkpts,[0:n-1])'.*[0.5:1:n-0.5].*pkw') * gfv # doing a legendre transform
  end
  mu == 0 && (gfv[2:2:end] = 0) # for centred transforms
  gfv /= gfv[1]*ds # to preserve probability density
  gfv
end

# Returns kernel matrix of a Fourier Gaussian distribution
spectralgauscoefs(n::Integer,periodic::Bool=false,dom::Array{Float64,2}=defdom(periodic),
                              sigma::Float64=defaultrelsigmasize*domsize(dom)[1],mu::Float64=0.) =
  periodic ? (fouriergauscoefs(n,dom,sigma)) : (leggauscoefs(n,dom,sigma))
spectralgauskernel(n::Integer,periodic::Bool=false,sratio::Float64=0.001) =
    spectralconv(spectralgauscoefs(n,periodic,defdom(periodic),domsize(defdom(periodic))[1]*sratio,0.),periodic,defdom(periodic),n)
spectralgauskernel(n::Integer,periodic::Bool=false,
                  dom::Array{Float64,2} = defdom(periodic),
                  sigma::Float64=defaultrelsigmasize*domsize(dom)[1]) =
  spectralconv(spectralgauscoefs(n,periodic,dom,sigma),periodic,dom,n)

legdeltacoefs(n::Integer,dom::Array{Float64,2}=defdom(false),ctr=mean(dom)) =
  DM.legp(ctr,[0:n-1],dom)[:].*[0.5:1:n-0.5]
fourierdeltacoefs(n::Integer,dom::Array{Float64,2}=defdom(true),ctr=dom[1]) =
  DM.fouriergauscoefs(n,dom,0.,ctr)
spectraldeltacoefs(n::Integer,periodic::Bool=false,dom::Array{Float64,2}=defdom(periodic),
                  ctr=(periodic ? dom[1] : mean(dom))) =
  periodic ? (fourierdeltacoefs(n,dom,ctr)) : (legdeltacoefs(n,dom,ctr))
