cd log
touch log-p-dummy
ls | grep log-p | xargs rm
cd ..
PIS="C2 D2 L2 M1 W3 W3 X3 X3 Y1 Y1"
for PI in $PIS; do
  bat log-p$PI julia -p 1 lri.jl $PI 40000 40000 const
done
