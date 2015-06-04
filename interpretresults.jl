Jpath = "/Users/johnwormell/Dropbox/Julia"
include("$Jpath/DiscreteMaps/DiscreteMaps.jl")
using HDF5, JLD, Dates, Gadfly, DiscreteMaps

searchdir(path,key) = filter(x->contains(x,key), readdir(path))
searchdirh5(path,key) = filter(x->contains(x,".h5"),searchdir(path,key))

function getresults(PInitial,Jpathx="",extrastuff="")
  path = "$(Jpath)$(Jpathx)"
  println("Looking in ",path)
  files = searchdirh5(path,"$(PInitial)$(extrastuff)")
  FL = length(files)
  println(FL, " files found")
  #  FL = 0 && (return (Nothing, Nothing))
  if FL > 0
    for i = 1:FL
      f = files[i]
      L = load("$(path)/$f")
      er = deepcopy(L["epsv"])
      eA = deepcopy(L["eA"])
      vA = deepcopy(L["vA"])
      if i == 1
        AN = size(eA,1)
        erv = deepcopy(er)

        eAv = deepcopy(eA)
        vAv = deepcopy(vA)
        println(size(erv),size(er),size(eAv),size(eA))

        println("AN: ",AN)
        println("$f contains ",length(er)," entries")
      elseif size(eA,1) == AN
        println("$f contains ",length(er)," entries")
        erv = deepcopy([erv, er])
        eAv = deepcopy([eAv  eA])
        vAv = deepcopy([vAv vA])
      end
      #    println(eAv[1:2,1:2])
      #    println(size(eAv))
      if i == FL
        #         evperm = sortperm(erv)[1:50]
        #        eAv = eAv[1,evperm]
        #                vAv = vAv[1,evperm]

        #                epsv = erv[evperm]


        #      println(size(erv),size(eAv))
        (erv,eAv,vAv) = DiscreteMaps.chopeps(erv,eAv,vAv,epsmax=0.0001)
        rss, pval, zeroval, lrtval = DiscreteMaps.checklinearresponse(erv,eAv,vAv,epsmax=0.0001)
        println("chi-square values for observables:",round(rss,5))
        println("P values for observables: ",round(pval,5))
        println("Max sd: ",sqrt(maximum(vAv[3,:])))
        return Gadfly.plot(x=erv,y=eAv[3,:],Geom.smooth,Geom.point),Gadfly.plot(x=1:length(pval),y=sort(pval))
      end
    end
  end
end

op, pp = getresults("L1","/results","")

op


