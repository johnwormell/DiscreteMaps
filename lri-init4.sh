ls | grep log- | xargs rm
bat log/log-C1 julia -p 4 lri.jl C2 40000000 40000
bat log/log-D1 julia -p 4 lri.jl D1 40000000 40000
bat log/log-L1 julia -p 4 lri.jl L1 40000000 40000
bat log/log-N1 julia -p 4 lri.jl M1 40000000 40000
