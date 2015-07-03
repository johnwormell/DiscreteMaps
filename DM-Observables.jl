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

function sino(n::Int64,m::Int64;phase::Float64=0.)
  A(x::Array{Float64}) = sin(2*pi*(n*x[1,:]+m*x[2,:])-phase)
end
function sino(n::Int64;phase::Float64=0.)
  A(x::Array{Float64}) = sin(2*pi*n*x-phase)
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
sin100A(phase::Float64=0.) = [sino(i;phase=phase) for i = 1:100]
sin100p30A = [sino(i;phase=phase) for phase=([0:14]/15 * pi), i = 1:100]
logiA = [ido(),sino(1),sino(100),gaussian(0.3,0.03),gaussian(0.3,0.10)]
logiA2 = [ido(),sino(1),sino(100),gaussian(0.4,0.0003),gaussian(0.4,0.003)]
coupA = [georg(100),jeroen(100)]
coupA2 = [gewgaw(100),jinx(100)]


Lhdefaultcpts = [0.95
 0.18050000000000013
 0.5620950500000004
 0.9353479781088903
 0.22979412423470444
 0.6725573818672567
 0.8368510098598475
 0.5188193091943236
 0.948654167685504
 0.1850958637100249
 0.5731744628003662
 0.9296528923767357
 0.2485138898747596
 0.7096679983734867
 0.7829494557406111
 0.6457705008851494
 0.869253652072432
 0.431876613638451
 0.9323649760764134
 0.23963000435728138
 0.6923883684022406
 0.8093495196733903
 0.5863509237758023
 0.9216653682596492
 0.27435360539972703
 0.756518077494812
 0.6999542084897893
 0.798069595127443
 0.6123871625501363
 0.9020026776369312
 0.3358966192564029
 0.8476663056283836
 0.49068693173670325
 0.949670413686188
 0.1816267724101116
 0.5648262542251891
 0.9340307156998984
 0.23414588375032924
 0.6814220377178049
 0.8249269680752137
 0.5488053685863065
 0.9409485367891881
 0.2111446740332013
 0.632937882510335
 0.8828445736959397
 0.39303412308798175
 0.9065215444704832
 0.32201288874906514
 0.8296182352684799
 0.5371369121182735]
Lhdefaultcmag = [0.5129891760425771
 0.2773927772077526
 0.1780134926840083
 0.2591299639225074
 0.14245987945869001
 0.09941190177042715
 0.08680901221424346
 0.0542549155101509
 0.14345987779648195
 0.07769048045232538
 0.05021943571353547
 0.06734188378629848
 0.03726658481243907
 0.02695602168403645
 0.021354181335058432
 0.014562022533437509
 0.013835029183363282
 0.00825868465841375
 0.011477729186155097
 0.0063317524063810525
 0.004501131075675547
 0.0037224207304112117
 0.002427693966583756
 0.0029967716314827854
 0.001674028648310165
 0.0012783267897511356
 0.0009155374338819563
 0.0007426839445294241
 0.0004934445997100272
 0.0005339163875520105
 0.0003054583361922263
 0.0002735183208502456
 0.0001682667567944215
 0.0006324775793625225
 0.00034213001922116216
 0.0002199460239800841
 0.00031335299887704854
 0.0001725308757476124
 0.0001213775413687152
 0.00010336809668929764
 6.57789360276959e-5
 0.00010800559402984843
 5.899911549560125e-5
 3.9819748701275964e-5
 3.961570055840278e-5
 2.3224656332232534e-5
 2.5758458664593332e-5
 1.4654509694748242e-5
 1.2599990543578592e-5
 7.96082006916964e-6]

peakedgeheight = 0.01
LhdefaultcptsA = [gaussian(Lhdefaultcpts[i],(Lhdefaultcmag[i]/peakedgeheight)^2) for i = 1:50]
sin100LhdefaultcptsA(phase::Float64=0.) = [sin100A(phase),LhdefaultcptsA]
