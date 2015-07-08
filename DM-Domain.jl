# Domain stuff

domsize(dom::Array{Float64,2}) = dom[:,2] - dom[:,1]
#domsize(M::Map) = domsize(M.dom)

restrictto(x::Array{Float64,1},xmin,xmax) = max(min(x,xmax),xmin)
function restrictto!(x::Array{Float64,1},xmax,xmin)
  x[:] = restrictto(x)
  nothing
end

function indomain(x::Array{Float64},dom::Array{Float64,2})
  return minimum(dom[:,1] .<= x .<=dom[:,2],1)
end
function checkindomain(x::Array{Float64},dom::Array{Float64,2})
  return x[:,vec(indomain(x,dom))]
end

# Distance from edge of (hyper-)rectangular domain
function domainedgedist(x::Array{Float64},dom::Array{Float64,2})
  return broadcast(min,abs(x.-dom[:,1]),abs(x.-dom[:,2]))
end

# Distance from edge of (hyper-)rectangular domain for x pointing in oriented directions
function domainedgedistoriented(x::Array{Float64},sgn::Array{Int8},dom::Array{Float64,2})
  dommatrix = (sgn .== -1) .* dom[:,1] + (sgn .== 1) .* dom[:,2]
  return abs(x .- dommatrix)
end

# Distance from nearest (hyper-rectangular) "barrier" for x pointing in oriented directions
function barrierdistoriented(x::Array{Float64},sgn::Array{Int8},dom::Array{Float64})
  distmatrix = fill(convert(Float64,-Inf),size(x))
  pmatrix = Array(Float64,length(x))
  for i = 1:size(dom,2)
    pmatrix[:] = (dom[:,i] .- x)[:].*sgn[:]
    distmatrix[:] = max(distmatrix[:],pmatrix)
  end
  return distmatrix
end

# Test function!
testfn(x::F64U,centre,width) = (abs(x.-centre) .< width) .* exp(1-(1 - ((x .-centre)./width).^2).^(-1))
normalisedtestfnspike(x::F64U) = (x .> 0) .* testfn(x,0,1) ./ sqrt(x)
normalisedtestfnspiketotalintegral = 1.526429236397948867946243218050187830844867513079036769761039753444952626979774
  # integral of testfn(x,0,1)/sqrt(x) on [0,1]
