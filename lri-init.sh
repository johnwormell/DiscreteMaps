ls | grep log- | xargs rm
bat log-C1 julia -p 2 lri.jl C1 40000000 40000
bat log-D1 julia -p 2 lri.jl D1 40000000 40000
bat log-L1 julia -p 2 lri.jl L1 40000000 40000
bat log-N1 julia -p 2 lri.jl N1 40000000 40000
bat log-X1 julia -p 7 lri.jl X1 40000000 40000
bat log-Y1 julia -p 7 lri.jl Y1 40000000 40000
bat log-W1 julia -p 7 lri.jl W1 40000000 40000