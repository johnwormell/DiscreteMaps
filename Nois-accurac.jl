include("DiscreteMaps.jl"); using DiscreteMaps; DM = DiscreteMaps

# Mp = DM.logistic(3.8)
Ns = [10:2:20,40:10:100,120:20:250,290:40:610]
L2err = Array(Float64,length(Ns))
timel = Array(Float64,length(Ns))
Ntop = 1000
# tic()
# DM.spectralacim(Mp,sigma=0.0,Ntop,verbose=true)
# toc()
# tic()
# rhotp = DM.spectralacim(Mp,sigma=0.01,Ntop,verbose=true).coeffs
# toc()

Mp = DM.spdoubling(0.02)
rhtp = DM.spectralacimpf(Mp,sigma=0.05,Ntop,verbose=true).coeffs
for i = 1:length(Ns)
  N = Ns[i]
  tic()
  rh = DM.spectralacimpf(Mp,sigma=0.05,N).coeffs
  coeffserr = (rhtp - [rh,zeros(Ntop-N)])
  L2err[i] = dot(coeffserr,DM.spectralinnerprodm(Ntop,Mp.periodic,Mp.dom)*coeffserr) |> sqrt
  timel[i] = toq()
  (log10(N), log10(L2err[i]),log10(timel[i]))|>println
end

[Ns log10(L2err)]
using PyPlot
PyPlot.plot(Ns,log2(L2err))
PyPlot.xlim([0,100])
#plot(x=Ns,y=L2err,Scale.x_log10,Scale.y_log10)
#plot(x=Ns,y=timel,Scale.x_log10,Scale.y_log10)

