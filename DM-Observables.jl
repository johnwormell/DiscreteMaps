export georg, jeroen, gj, sino, coso, ido, trigA, logiA, coupA

function georg(n::Int64)
  A(x::Array{Float64}) = sin(sum(sin(n*2*pi*x),1)/size(x,1))
end
function jeroen(n::Int64)
  A(x::Array{Float64}) = sum(sin(n*2*pi*x),1)/sqrt(size(x,1))
end
function gewgaw(n::Int64)
  A(x::Array{Float64}) = sin(sum(sin(n*2*pi*x),1)/size(x,1))
end
function jinx(n::Int64)
  A(x::Array{Float64}) = sum(sin(n*2*pi*x),1)/(size(x,1))
end
function gewgawaux(n::Int64)
  A(x,::Array{Float64}) = sin(sum(sin(n*2*pi*x),1)/sqrt(size(x,1)))
end

function gj(n::Int64,t::Int64,coup=true)
  t==1 && return georg(n), 1/2
  t==2 && return jeroen(n), 0
  t==3 && return gewgaw(n), 0
  t==4 && return jinx(n), 0
  t==5 && (coup==true ? (return gewgawaux(n), 0) : (return gewgaw(n), 1/2))
  t==6 && return jinx(n), 1/2
end

function sino(n::Int64,m::Int64)
  A(x::Array{Float64}) = sin(2*pi*(n*x[1,:]+m*x[2,:]))
end
function sino(n::Int64)
  A(x::Array{Float64}) = sin(2*pi*n*x)
end
function coso(n::Int64,m::Int64)
  A(x::Array{Float64}) = cos(2*pi*(n*x[1,:]+m*x[2,:]))
end
function coso(n::Int64)
  A(x::Array{Float64}) = cos(2*pi*n*x)
end
function ido(dim::Int64=1)
  A(x::Array{Float64}) = x[dim,:]
end

function gaussian(mu,sd)
  A(x::Array{Float64}) = exp(-(x-mu).^2/(2*sd^2))/sd/sqrt(2pi)
end

trigA = [ido(),sino(1),sino(2),sino(3),sino(10),sino(100),sino(1000),coso(1),coso(2),coso(3),coso(10),coso(100),coso(1000)]
trig2A = [ido(1),ido(2),sino(1,0),sino(0,1),sino(3,2),sino(10,0),sino(100,0),sino(1000,0),coso(0,1),coso(3,2),coso(10,0),coso(100,0),coso(1000,0)]
logiA = [ido(),sino(1),sino(100),gaussian(0.3,0.03),gaussian(0.3,0.10)]
coupA = [georg(100),jeroen(100)]
coupA2 = [gewgaw(100),jinx(100)]
