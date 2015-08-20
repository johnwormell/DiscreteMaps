function spectralacimpf(M::IMap, # map whose acim we are finding
                      N::Integer=100; # number of points to take spectra at
                      verbose=false, # print output about accuracy etc
                      returntransfmat=false, # if no critical points, additionally return the transfer operator matrix
                      sigma::Float64 = 0.)
  noiseexists = (sigma > 0.)

  the_spectralpts = spectralpts(N,M.periodic,M.dom)

  U = spectralf(M.f(the_spectralpts,M.params),[0:N-1], M.periodic,M.dom) # Koopman operator
  U = spectraltransf(N,M.periodic)*U
  chopm!(U)
  L = inv(spectralinnerprodm(N)) * U' * spectralinnerprodm(N)
  if noiseexists
    spectralker = spectralgauskernel(N,M.periodic,M.dom,sigma) #|> chopm
    L = spectralker * L
  end
  (r,(U,S,V)) = findinvmeasure(L,verbose=verbose)

  mu = SpectralMeasure(r,M.dom,M.periodic)

  # Normalising
  munorm = totalint(mu)[1]
  mu = mu/munorm

  # checking everything's ok
  if verbose
    println("r: ",round(r[1:10]'/munorm,4))
    println("Lr: ",round((L*r)[1:10]'/munorm,4))
  end

  return returntransfmat ? (mu, L) : mu
end
