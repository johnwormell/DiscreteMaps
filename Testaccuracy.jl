using HDF5, JLD

L = JLD.load("/Users/johnwormell/Dropbox/Julia/Specacim/plotoutput.h5")

using Gadfly
lsp = linspace(0.05,0.95,length(L["bgnorm"]))
Gadfly.plot(layer(x=lsp,y=L["bgnorm"]/5,Geom.line,Theme(default_color=color("orange"))),
            layer(x=lsp,y=L["xhcts"],Geom.bar),Scale.y_continuous(maxvalue=4))

difference = (L["xhcts"]./L["bgnorm"]*5 - 1)
goodvals = abs(difference) .< 0.1
sqrt(L["xhcts"][goodvals])
mean(L["xhcts"] / (1/ 10^7 * 10^3 / 0.9))
Gadfly.plot(layer(x = lsp[goodvals],y=difference[goodvals].*sqrt(L["xhcts"][goodvals]/ (1/ 10^7 * 10^3 / 0.9)),Geom.histogram))#,Scale.y_continuous(maxvalue=0.1))


sum(difference[goodvals].*sqrt(L["xhcts"][goodvals]/ (1/ 10^7 * 10^3 / 0.9)).^2)

length(goodvals)

