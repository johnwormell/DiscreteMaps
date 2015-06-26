cd log
touch log-s-dummy
ls | grep log-s | xargs rm
cd ..
PIS="C2 2
D2 2
L2 2
M1 2
W3 7
X3 7
Y1 7
"
for PI in $PIS; do
  set -- $PI
  bat log-s$1 julia -p $2 lri.jl $1 40000000 40000
done
