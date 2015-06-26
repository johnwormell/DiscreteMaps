cd log
touch log-s-dummy
ls | grep log-s | xargs rm
cd ..
PIS="C2 D2 L2 M1 N1 W3 W3 W3 X3 X3 X3 Y1 Y1 Y1
"
for PI in $PIS; do
  bat log-s$PI julia -p 2 lri.jl $PI 40000000 40000
done
