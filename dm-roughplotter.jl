include("DiscreteMaps.jl")
using DiscreteMaps, Distributions, Gadfly
using HDF5, JLD
#W = load("Dropbox/Julia/results/rs-W1-2015-06-03T10-00-00--15-5.h5")
#X = load("Dropbox/Julia/results/rs-X1-2015-06-03T10-00-00--15-5.h5")
#Y = load("Dropbox/Julia/results/rs-Y1-2015-06-03T10-00-00--15-5.h5")
C = load("/Users/johnwormell/Dropbox/Julia/DiscreteMaps/results/rp-L1-2015-06-05T13-00-00--18-26.h5")
DiscreteMaps.checklinearresponse(C["epsv"],C["eA"],C["vA"])
maximum(sqrt(C["vA"][2,:]))

#using Interact
i=1
Gadfly.plot(
  x=C["epsv"],y=C["eA"][i,:],Geom.point,Theme(default_color=color("blue")),Geom.smooth(method=:lm)
#  layer(x=C["epsv"],y=C["eA"][2,:],Geom.point,Theme(default_color=color("red"))),
#        layer(x=Y["epsv"],y=Y["eA"][2,:],Geom.point,Theme(default_color=color("green")))
  )

g
xh, m = DiscreteMaps.obsiterate(DiscreteMaps.doubling1(N=10^5),)

#Gadfly.plot(x=[1:100],y=xh[10000:11100])

