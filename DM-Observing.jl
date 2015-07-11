export observe, autocovariance, observevar

function observe(A::Function,x_history::Array{Float64})
  # Estimates an observable
  return mean(A(x_history))::Float64#::Array{Float64} #, var(A(x_history))/length(x_history)::Array{Float64}
end

function autocovariance(A::Function,x_history::Array{Float64},sumN = 60)
  NH = length(x_history)
  A1 = A(x_history)[:]
  A2 = A(x_history)[:]
  varA = Array(Float64,sumN+1)
  for i = 0:sumN
    varA[i+1] = cov(A1,A2) #(sum(A1.*A2)-sum(A1)*sum(A2)/(N-i))/(N-i-1)
    #        flucterms[i+1] = mean(flucintegrand)
    #        flucvar[i+1] = var(flucintegrand)/(N-i)
    shift!(A1)
    pop!(A2)
  end
  return varA
end
function observevar(A::Function,x_history::Array{Float64},sumN = 60)
  # Estimates the random variance of a measured observable
  varA = autocovariance(A,x_history,sumN)
#  println(cov(A1,sort(A1)))
#  println(floor(varA))
  varA *= 2
  varA[1] /= 2
  cesaroseries = cumsum(varA)
  return sum(cesaroseries)/(sumN+1)
end
