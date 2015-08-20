include("DiscreteMaps.jl"); using DiscreteMaps; DM = DiscreteMaps

epsv = [-1e-5:2e-7:1e-5]
Neps = length(epsv)
Nnonoise = 400

idor = Array(Float64,Neps)
for i = 1:Neps
  idor[i] = DM.totalint(DM.spectralacim(DM.logistic(3.8*(1+epsv[i])),Nnonoise,uselogisticcofn=true),identity)[1]
end
sigmav = [9:-2:1]/5000
Nsigma = length(sigmav)

idos = Array(Float64,Nsigma)
idofds = Array(Float64,Nsigma)
Mp2 = DM.logistic(3.8); Mp2.dom = [0 1.]

idoconst = DM.spectralinnerprodm(100,Mp2.periodic,Mp2.dom)[2,2]
for i = 1:Nsigma
  (mu,L) = DM.spectralacimpf(Mp2,int(3/sigmav[i]),sigma=sigmav[i],returntransfmat=true,verbose=true)
  idos[i] = mu.coeffs[2]/idoconst
  idofds[i] = DM.flucdis(mu,L,[0.,1,zeros(length(mu.coeffs)-2)]).coeffs[2]/idoconst
end

sigmav
