cd log
ls | grep log-p | xargs rm
cd ..
bat log-pC2 julia -p 2 lri.jl C2 40000 40000 const
bat log-pD2 julia -p 2 lri.jl D2 40000 40000 const
bat log-pL2 julia -p 2 lri.jl L2 40000 40000 const
bat log-pM1 julia -p 2 lri.jl M1 40000 40000 const
bat log-pW3 julia -p 7 lri.jl W3 40000 40000 const
bat log-pX3 julia -p 7 lri.jl X3 40000 40000 const
bat log-pY1 julia -p 7 lri.jl Y1 40000 40000 const