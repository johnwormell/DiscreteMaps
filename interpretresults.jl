Jpath = "/Users/johnwormell/Dropbox/Julia"
include("$Jpath/DiscreteMaps/DiscreteMaps.jl")
using HDF5, JLD, Dates, Gadfly, DiscreteMaps

#end

function plotresults(PInitial,Jpathx="",extrastuff=["rs"])
  erv, eAv, vAv = DiscreteMaps.synthesiseresults(PInitial,Jpathx,extrastuff)
        (erv,eAv,vAv) = DiscreteMaps.chopeps(erv,eAv,vAv,epsmax=0.0001)
        rss, pval, zeroval, lrtval = DiscreteMaps.checklinearresponse(erv,eAv,vAv,epsmax=0.0001)
        println("chi-square values for observables:",round(rss,5))
        println("P values for observables: ",round(pval,5))
        println("Max sd: ",sqrt(maximum(vAv[3,:])))
        return Gadfly.plot(x=erv,y=eAv[3,:],Geom.smooth,Geom.point),Gadfly.plot(x=1:length(pval),y=sort(pval))
end

op, pp = plotresults("L1","/Users/johnwormell/Dropbox/Julia/DiscreteMaps/results",["rs","2015"])

typeof(op)


