include("DiscreteMaps.jl")
using DiscreteMaps, Gadfly
println("Starting")
M = DiscreteMaps.logistic()
CO = DiscreteMaps.criticalorbit(M)
Sp = DiscreteMaps.Spikes(CO,M.dom,[1.])

norminv(Mat) = 1/minimum(svdvals(Mat))
crit = M.crit(M.params)[1]
N = 1000
L = Array(Float64,N,N)
chebypts = DiscreteMaps.chebypts(N,M.dom)
chebypts[N] -= 2eps(1.)
chebypts[1] += 2eps(1.)
V = DiscreteMaps.chebytn(chebypts,[1:N],crit,M.dom)

for i = 1:N
  L[:,i] = DiscreteMaps.transfer(DiscreteMaps.chebytn,M, (i,crit,M.dom))(chebypts)
end # only T_1 ... T_N cause it's got subtraction happening
#L[N,:] = L[N-1,:]
L
p = (DiscreteMaps.spikefn(chebypts,Sp) - DiscreteMaps.transfer(DiscreteMaps.spikefn,M,([Sp]))(chebypts))
#p[N] = p[N-1]
rcheck = (DiscreteMaps.shortchebytransf(N) * L - I) \ (DiscreteMaps.shortchebytransf(N) * p)
pc = DiscreteMaps.spikefn(chebypts,Sp)
#(DiscreteMaps.transfer(DiscreteMaps.spikefn,M,([Sp]))(chebypts)./pc)[900:1000]
norminv((DiscreteMaps.shortchebytransf(N) * L - I))
scl(x) = (x) .* (100 .> abs(x))
Gadfly.plot(layer(x=chebypts,y=(V*rcheck-(V*rcheck)[200])/10/N^2,Geom.line,Theme(default_color=color("orange"))),
            layer(x = chebypts,y=(pc[1:990]),Geom.line))
#PyPlot.plot(chebypts,V*rcheck)

using HDF5, JLD
Ld = load("/Users/johnwormell/Dropbox/Julia/DiscreteMaps/results/rd-L1-2015-06-04T18-56-26.h5")


println("Done")
