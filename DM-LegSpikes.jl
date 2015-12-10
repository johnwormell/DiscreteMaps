# Multiplying
function legnakedspikecoeffs(Sp::Spikes,n::Integer,dom::Array{Float64,2}=Sp.dom;
                             normspikes=true) # false = multiply by mag/mag0
  Nc = Sp.CO.Nc
  Npts = Sp.CO.Npts
  domsizemult = sqrt(2/domsize(dom)[1])
  coefm = Array(Float64,n,Npts,Nc)
  for i = 1:Nc
    for k = 1:Npts
      csgn = Sp.CO.sgn[k,i]
      cnorm = csgn * legnormcoords(Sp.CO.pts[k,i],dom)
      csgn2 = cnorm >= 0 ? 1 : -1
      theta = acos(cnorm*csgn2)
      coefm[:,k,i] = domsizemult * sqrt(2) *
        (csgn2 == 1 ? sin : cos)(theta*collect(1/2:1:n-1/2))
      normspikes || (coefm[:,k,i] .*= Sp.mag0[i] .*Sp.CO.mag[k,i])
      coefm[:,k,i] .*= (csgn2 * csgn).^ collect(0:n-1)

      #       coefm[1,k,i] = Sp.mag0[i]*Sp.CO.mag[k,i] * domsizemult * 2sqrt(1-cnorm)
      #       coefm[2,k,i] = csgn * coefm[1,k,i] * (cnorm + (1-cnorm)/3)
      #       for j = 1:n-2
      #         coefm[j+2,k,i] =  (csgn*2cnorm*(2j+1) * coefm[j+1,k,i] - (2j-1) * coefm[j,k,i])/(2j+3)
      #       end
    end
  end
  return coefm# .*[1:2:2n-1]
end

function legclothingspikecoeffs(Sp::Spikes,n::Integer,dom::Array{Float64,2}=Sp.dom;
                                normspikes=true)
  Nc = Sp.CO.Nc
  Npts = Sp.CO.Npts
  if Sp.widths == nothing
    return zeros(n,Nc,Npts)
  else
    coefm = Array(Float64,n,Npts,Nc)
    L = legtransf(n)
    Lpts = legpts(n,dom)
    for i = 1:Nc
      for k = 1:Npts
        coefm[:,k,i] = - L * spikefn(Lpts,Sp,k,nothing,i,oneminustfn=true)
        normspikes && (coefm[:,k,i] /= Sp.mag0[i] * Sp.CO.mag[k,i])
      end
    end
    return coefm
  end
end

function legspikecoeffs(Sp::Spikes,n::Integer,dom::Array{Float64,2}=Sp.dom;
                        normspikes=true)
  coefm = legnakedspikecoeffs(Sp,n,dom,normspikes=normspikes)
  if Sp.widths != nothing
    coefm += legclothingspikecoeffs(Sp,n,dom,normspikes=normspikes)
  end
  coefm
end

function legspikeinnerprodm(mu::SumMeasure,normspikes=true) # L_infty x L_1 -> R
  n = mu.components[1].N
  [diag(leginnerprodm(n,mu.dom)) .* legspikecoeffs(mu.components[2],n,mu.dom,normspikes=true)[:,:] full(leginnerprodm(n,mu.dom))]
end

function legspikemult(coeffs::Array{Float64,1}, Sp::Spikes, n::Integer=length(coeffs), dom = Sp.dom; spikesonly=false)
  Sp.widths == nothing || error("legspikemult doesn't support spikes that have test functions on them")
  spn = Sp.CO.Nc * Sp.CO.Npts

  SL = legnakedspikecoeffs(Sp, n+1, dom, normspikes=true)[:,:] # reusing SL as spike Legendre coeffs matrix to save memory
  for i = 1:spn
    lc = copy(SL[:,i])
    lp = legp(Sp.CO.pts[i],collect(0:n),dom)

    mconstarray2 = zeros(n+1)
    SL[:,i] = coeffs[1] * mconstarray2
    mconstarray1 = [lc[2]/3-lp[2]*lc[1],
                    -0.5*(lc[3:end]./collect(5:2:2n+1) - lc[1:end-2]./collect(1:2:2n-3)),
                    0]
    coeffs[2] != 0 && (SL[:,i] += coeffs[2] * mconstarray1)
    for j = 2:maximum([find(coeffs.!=0),2])-1
      mconstarray = -mconstarray2*(j-1)/j
      mconstarray[1] = lc[j+1]/(2j+1) - lc[1] * lp[j+1]
      mconstarray[end-1] = 0

      for k = 1:n-2
        mconstarray[k+1] += (2j-1)/j * (k+1)/(2k+3) * mconstarray1[k+2]
      end
      for k = 1:n-1
        mconstarray[k+1] += (2j-1)/j * (k/(2k-1) * mconstarray1[k] -
                                          lp[j]/2 *(lc[k+2]/(2k+3) - lc[k]/(2k-1)))
      end

      coeffs[j+1] != 0 && (SL[:,i] += coeffs[j+1] * mconstarray)
      mconstarray2 = mconstarray1
      mconstarray1 = mconstarray
    end
  end

  SL /= domsize(dom)[1]
  if spikesonly
    SS = legapprox(vec(Sp.CO.pts),coeffs,dom)
    Sp.widths != nothing && error("Can't do spikes only if the spikes have test functions on them")
    return (SS, SL[1:n,:])
  end
  LL = legmult(coeffs, n)
  SS = legapprox(vec(Sp.CO.pts),coeffs,dom)

  if Sp.widths != nothing
    SL[1:n,:] += LL*legclothingspikecoeffs(Sp,n,dom,normspikes=true)[:,:]
  end

  multm = [diagm(SS) zeros(spn,n);
           SL[1:n,:] LL]
end

function legspikelinop(L::Array{Float64,2}, mu::SumMeasure)
  Sp = deepcopy(mu.components[2])
  Sm = deepcopy(mu.components[1])
  sn = L * [vec(Sp.mag0' .* Sp.CO.mag),Sm.coeffs]
  Sp.mag0 = ones(size(Sp.mag0));
  spn = Sp.CO.Nc * Sp.CO.Npts
  Sp.CO.mag[:] = sn[1:spn]
  Sm.coeffs = sn[spn+1:end]
  Sm + Sp
end

function legspikesum(coeffs::Array{Float64,1}, mu::SumMeasure)
  Sp = deepcopy(mu.components[2])
  Sm = deepcopy(mu.components[1])
  Sp.mag0 = ones(size(Sp.mag0));
  spn = Sp.CO.Nc * Sp.CO.Npts
  Sp.CO.mag[:] = coeffs[1:spn]
  Sm.coeffs = coeffs[spn+1:end]
  Sm + Sp
end

function legspikegetcoeffs(mu::SumMeasure)
  Sp = mu.components[2]
  Sm = mu.components[1]
  [(Sp.CO.mag .* Sp.mag0')[:],Sm.coeffs]
end

function legspikecorrel(L::Array{Float64,2},mu::SumMeasure)
  Sp = mu.components[2]
  if Sp.widths == nothing
    Spint = sqrt(domsize(mu.dom)[1]/2) * 2sqrt(1 - Sp.CO.sgn .* legnormcoords(Sp.CO.pts,mu.dom))
  else
    Spint = normalisedtestfnspiketotalintegral .* sqrt(Sp.widths) |> vec # assuming mag0 = mag = 1 for the transfer operator
  end

  n = mu.components[1].N
  spn = Sp.CO.Nc * Sp.CO.Npts

  invs = inv(I-L[[1:spn,spn+2:spn+n],[collect(1:spn),collect(spn+2:spn+n)]])
  invs = [invs[1:spn,:];
          -Spint' * invs[1:spn,:] / legtotalint(n,mu.dom)[1];
          invs[spn+1:end,:]]
  invs = [invs[:,1:spn] zeros(size(invs,1)) invs[:,spn+1:end]]
  #   inv = [invs[1:spn,1:spn] zeros(spn) invs[1:spn,spn+2:end];
  #           zeros(1,spn+n);
  #           invs[spn+2:end,1:spn] zeros(n-1) invs[spn+2:end,spn+2:end]]
  #   Spint = normalisedtestfnspiketotalintegral * Sp.mag0 .* Sp.CO.mag .* sqrt(Sp.widths) |> vec
  #   invf[spn+1,1:spn] = -Spint' * invs[:,1:spn] / legtotalint(n,mu.dom)[1]
  return invs
end

function legspiketoleg(Sp::Spikes,n)
  [legnakedspikecoeffs(Sp,n,Sp.dom+100eps(1.)*[-1 1.],normspikes=true)[:,:] eye(n)]
end# function legspikemult(coeffs::Array{Float64,1},Sp::Spikes,n::Integer)
