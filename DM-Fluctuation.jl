#export fluctuation, flucsum
function fluctuation(M::Map,A::Function,x_history::Array{Float64},divrX::Function,flucN=60)
  N = length(x_history)
  divrXxh = divrX(x_history)[:]
  Axh = A(x_history)[:]

  flucterms = Array(Float64,flucN+1)
  flucvar = Array(Float64,flucN+1)
  for i = (0:flucN)
    flucintegrand = - divrXxh.*Axh
    flucterms[i+1] = mean(flucintegrand)
    flucvar[i+1] = var(flucintegrand)/(N-i)
    shift!(Axh)
    pop!(divrXxh)
  end
  return flucterms, flucvar
end

function fluctuation(M::Map,A::Array{Function},x_history::Array{Float64},divrX::Function,flucN=60)
  AN = length(A)
  flucterms = Array(Float64,flucN+1,AN)
  flucvar = Array(Float64,flucN+1,AN)
  for i = 1:AN
    Ai = A[i]
    ft, fv = fluctuation(M,Ai,x_history,divrX,flucN)
    flucterms[:,i] = ft
    flucvar[:,i] = fv
  end
  return flucterms, flucvar
end

# This method is crap for runtime don't use it it seems
## function fluctuation(M::Map,A::Array{Function},x_history::Array{Float64},divrX::Function)
##     flucN = 60
##     AN = length(A)
##     flucterms = Array(Float64,flucN+1,AN)
##     flucvar = Array(Float64,flucN+1,AN)

##     N = length(x_history)
##     divrXxh = divrX(x_history)[:]
##     Axh = Array(Float64,N,AN)
##     for j = 1:AN
##         Axh[:,j] = A[j](x_history)[:]
##     end

##     flucintegrand = Array(Float64,N,AN)
##     for i = (0:flucN)
##         for j = 1:AN
##             flucintegrand[:,j] = - divrXxh.*Axh[:,j]
##         end
##         flucterms[i+1,:] = mean(flucintegrand,1)*N/(N-i)
##         flucvar[i+1,:] = var(flucintegrand,1)*(N-1)/(N-i-1)/(N-i)

##         pop!(divrXxh)
##         unshift!(divrXxh,0)
##     end
##     return flucterms, flucvar
## end

# F-D summation

function flucsum(flucterms::Array,flucvar::Array,limsuplength=10,output=true)
  NA = size(flucterms,2)
  flucN = size(flucterms,1) - 1

  p = 1e-7
  limsupran = (flucN-limsuplength+2:flucN+1)
  ftlimsupran = flucterms[limsupran,:]
  # Chi-square test on the tail of flucterms to check if there is a pole at lambda = 1
  polenormalised = (ftlimsupran - repmat(mean(ftlimsupran,1),limsuplength,1))./sqrt(flucvar[limsupran])

  if minimum(Distributions.ccdf(Distributions.Chisq(limsuplength-1),sum(polenormalised.^2,1))) > p/NA
    termsconverge = true
  else
    termsconverge = false
  end

  # Chi-square test on the tail of flucterms to check if it is non-zero
  divnormalised = ftlimsupran./sqrt(flucvar[limsupran])
  if true# Distributions.ccdf(Distributions.Chisq(limsuplength), sum(divnormalised.^2)) < p
    return sum(flucterms,1),sum(flucvar,1)
  else
    termsconverge && output && print("There might be a pole at lambda = 1")
    cesaroseries = cumsum(flucterms)
    cesarovariation = cumsum(flucvar)
    return mean(cesaroseries),mean(cesarovariation)/(flucN+1),termsconverge
  end
  # TODO:
  ## make first chi-square test exact
  ## make it return Float and not Array for 1-dimensional inputs
  ## estimate expectation and variance for background ind. variables -> know when terms converge
  ## proper analytic continuation
end
