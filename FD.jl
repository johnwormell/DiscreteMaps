using Gadfly
include("DiscreteMaps.jl");using DiscreteMaps; DM = DiscreteMaps
# Mp = DM.logistic()
# (mu,Lhat) = DM.spectralacim(Mp,200,returntransfmat=true)
# mu = mu.components[1]
# Lhat = Lhat[2:end,2:end]
# DMmu.components[1].coeffs
function gaussianfouriercoefs(n::Integer,dom::Array{Float64,2}=DM.defdom(true),sigma::Float64=0.001*DM.domsize(dom)[1],mu::Float64=0.)
  ds = DM.domsize(dom)[1]
  omega0 = mu*2pi/ds
  gfv = Array(Float64,n)
  gfv[1] = 1. / ds
  for i = 2:n
    i2 = div(i,2)
    gfv[i] = (rem(i,2)==0 ? cos(-i2*omega0) : sin(-i2*omega0)) *
      exp(-(2pi/ds * sigma * i2)^2 / 2) * (2 / ds)
  end
  gfv
end
function fouriergauskernelm(n::Integer,sratio::Float64=0.001) = #dom = DM.defdom(true),sigma::Float64=0.001*DM.domsize(dom)[1])
  DM.fourierconv(gaussianfouriercoefs(n,[0. 2pi],2pi*sratio,0.),[0. 2pi],n)
function fouriergauskernelm(n::Integer,dom = DM.defdom(true),sigma::Float64=0.001*DM.domsize(dom)[1])
  DM.fourierconv(gaussianfouriercoefs(n,dom,2pi*sigma,0.),[0. 2pi],n)

periodic=true;n = 100
fpts = DM.fourierpts(n,[0. 2pi]);
  plot(x=fpts,y=DM.fourierapprox(fpts, gaussianfsc(n,[0. 2pi]),[0. 2pi]),Geom.line)


DM.spectralapprox(DM.spectralpts(n,true),[1. ])

a = 0.5.^[1:30]; b = 0.8.^[1:30]
[DM.fourierconv(a)*b DM.fourierconv(b)*a]
    Mp = DM.spdoubling(0.1,1)
Mp.periodic = false
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

