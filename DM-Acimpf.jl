function spectralacimpf(M::Map, # map whose acim we are finding
                        N::Integer=100; # number of points to take spectra at
                        verbose=false, # print output about accuracy etc
                        returntransfmat=false, # if no critical points, additionally return the transfer operator matrix
                        sigma::Float64 = 0., # width of noise
                        fastinv = false) # options: false, :randn, :ortho
  noiseexists = (sigma > 0.)

  the_spectralpts = spectralpts(N,M.periodic,M.dom)

  Pf = spectralf(M.f(the_spectralpts,M.params),[0:N-1], M.periodic,M.dom) # Koopman operator
  Pf = spectraltransf(N,M.periodic)*Pf
  chopm!(Pf)
  L = inv(spectralinnerprodm(N,M.periodic,M.dom)) * Pf' * spectralinnerprodm(N,M.periodic,M.dom)
  if noiseexists
    spectralker = spectralgauskernel(N,M.periodic,M.dom,sigma) #|> chopm
    L = spectralker * L
  end
  if fastinv == false
    (r,(U,S,V)) = findinvmeasure(L,verbose=verbose)
  elseif fastinv == :rand
    # uses random vector and gap in svd values of (L-I) to look into the kernel of (L-I).
    # Assumes (L-I) non-singular, which it really should be. In an orthonormal basis you
    # can make it so by using chopm
    rn = rand(N)
    r = (L-I)\rn
    verbose && println("Smallest singular value is approximately ",norm(rn)/norm(r))
  elseif fastinv == :ortho
    # uses separation of constant function and its L2 orthogonal complement,
    # i.e. that the top row of L-I is empty
    norm(chopm((L-I)[1,:]|>vec)) > 0 && warn("Non-zero entries in the top column of (L-I)...")
    r = [1,-(L-I)[2:end,2:end]\(L-I)[2:end,1]]
  end


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
