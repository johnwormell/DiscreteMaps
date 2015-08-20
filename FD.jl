using Gadfly
include("DiscreteMaps.jl");using DiscreteMaps; DM = DiscreteMaps

# Mp = DM.logistic()
# (mu,Lhat) = DM.spectralacim(Mp,200,returntransfmat=true)
# mu = mu.components[1]
# Lhat = Lhat[2:end,2:end]
# DMmu.components[1].coeffs


Mp = DM.spdoubling(0.1,3)
#Mp.periodic = false
(mu, Lhat) = DM.spectralacim(Mp,100,returntransfmat=true)
#mu.periodic = false

X(x) = -0.02 * sin(2*pi*x) #x .* (1 - x)
function dq(x)
  fx = mod(2*x+2*0.1*sin(2*pi*x),1)
  return fx + X(fx) #0.002*cos(2*pi*2*fx)
end
xh = DM.iterate(DM.desque(dq),10^6)
#xh = DM.iterate(DM.spdoubling(0.02/2,1),10^6)
dmu = DM.flucdis(mu,Lhat,X)
(pgd,dmud) = DM.plotmeasure(dmu)
mud = DM.plotmeasure(mu)[2]
plot(layer(x=pgd,y=(mud+dmud)*10^6/200,Geom.line,Theme(default_color=color("orange"))),
     layer(x=pgd,y=mud*10^6/200,Geom.line,Theme(default_color=color("red"))),
     layer(x=xh,Geom.histogram(bincount=200))
     )

