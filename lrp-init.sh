cd log
touch log-p-dummy
ls | grep log-p | xargs rm
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
  bat log-p$1 julia -p $2 lri.jl $1 40000 40000 const
done
