ls | grep log- | xargs rm
bat log-C1-2 julia -p 2 lri.jl C1 40000000 40000
bat log-M1 julia -p 4 lri.jl L1 40000000 40000
bat log-X2 julia -p 14 lri.jl X2 40000000 40000
bat log-W2 julia -p 14 lri.jl W2 40000000 40000