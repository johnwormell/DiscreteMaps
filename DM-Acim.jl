
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
    pts[1,i] = M.f(crit[i],M.params)[1]
    mag[1,i] = sqrt(2/abs(critxx[i]))
    sgn[1,i] = sign(critxx[i])

    for p = 1:(Npts-1)
      pts[p+1,i] = M.f(pts[p,i],M.params)[1]
      dfp = (M.df(pts[p,i],M.params))[1]
      mag[p+1,i] = mag[p,i] ./ sqrt(abs(dfp))
      sgn[p+1,i] = sgn[p,i] * sign(dfp)
    end
  end
  # periodic error?
  return CriticalOrbit(crit,pts,mag,sgn)
end

# Computes a function containing all the spikes from the map.
# The acim should be continuous after these spikes are removed.
function spikefn(x::Array{Float64,1},Sp::Spikes)
  Nx = length(x)
  Nc = Sp.CO.Nc
  xd = [x,Sp.CO.crit]'
  rawspikes = Array(Float64,Sp.CO.Npts,Nx+Nc)
  prel = zeros(Nx+Nc)
  for i = 1:Nc
    rawspikes[:,:] = (Sp.mag0[i] *
                        # size of initial bump
                        (sign(xd.-Sp.CO.pts[:,i]) .== Sp.CO.sgn[:,i]) .*
                      # is the x value on the right side of the spike?
                      abs(Sp.CO.mag[:,i]) # relative total size of spike vs bump
                      .* abs(xd.-Sp.CO.pts[:,i]).^(-0.5) # size of spike at x value
                      )
    nonedgecpts = Sp.CO.pts[(Sp.CO.pts .!= Sp.dom[:,1]) & (Sp.CO.pts .!= Sp.dom[:,2])] # only works in 1D
    prel += (indomain(xd,Sp.dom) .*
             # is the x value in the domain of the map?
             sum(rawspikes .*
                 testfn(xd,Sp.CO.pts[:,i],
                        # test function centred around spikes to limit their size
                         min(Inf,#minimum(domainedgedist(nonedgecpts,Sp.dom)), #cbf
                            Sp.width[i])
                        # choice of width keeps them away from domain edges
             ),1)) |> vec

  end
  bumpwdth = minimum([diff(sort(Sp.CO.crit[:])),Sp.width])
  bumpsize = Sp.mag0 - prel[Nx+(1:Nc)]
  for i = 1:Nc
    prel += sum(testfn(xd,Sp.CO.crit[i],bumpwdth) .* bumpsize,1) |> vec # this will break for more than one critical point if they're too close to each other :/
  end
  return prel[1:Nx]

end
spikefn(x::Float64,Sp::Spikes) = spikefn([x],Sp)

# Chebyshev functions and transforms - 1D ONLY
inormcoords(x::F64U,dom::Array{Float64}=[-1.,1.]) = (dom[1] + dom[2])/2 + x * (dom[2] - dom[1])/2  #  - 3*eps(dom[2]-dom[1])
normcoords(x::F64U,dom::Array{Float64}=[-1.,1.]) = (x - (dom[1] + dom[2])/2)/(dom[2] - dom[1]) * 2

chebyt(x::F64U,k::I64U,dom::Array{Float64}=[-1.,1.]) =cos(acos(normcoords(x,dom)) * k')
chebytn(x::F64U,k::I64U,critp::Float64, dom::Array{Float64}=[-1.,1.]) = chebyt(x,k,dom) .- chebyt(critp,k,dom)
chebypts(n::Integer,dom::Array{Float64}=[-1.,1.]) = inormcoords(cos(linspace(-pi,0,n)),dom)
chebytransf(n::Integer) =
  diagm([1,2*ones(n-1)]) * cos([0:(n-1)]*linspace(-pi,0,n)') *
  diagm([1,2*ones(n-2),1]) / (2n-2)
shortchebytransf(n::Integer) = cos([1:(n)]*linspace(-pi,0,n)') .* [1,2*ones(n-2),1]' / (n-1)
# function chebyintegrate(v::Vector{Float64})
#   N = length(v)
#   v = [0,v[2:N]./[1:N-1],0] - [v[1:N-1]


# Transfer operator aka Frobenius-Perron operator
function transfer(r::Function,M::IMap,rargs=())
  function Lr(x::Array{Float64})
    Lrx = Array(Float64,size(x))
    for i = 1:length(x)
      Lrx[i] = sum(r(vec(M.g(x[i],M.params)),rargs...) ./
                   abs(vec(M.df(M.g(x[i],M.params),M.params))))
    end
#    Lrx[:] += r(M.g[i](x[:])) .* broadcast(M.df(M.g[i](r)))
    return Lrx
  end
  Lr(x::Float64) = Lr([x])[1]
  return Lr
end




