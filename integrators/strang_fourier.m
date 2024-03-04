function U = strang_fourier(U,M,gflux,T,m,ft,ift)
  tau = T/m;

  EE = exp(tau*M);
  U = gflux(tau/2,U);
  for jj = 1:(m-1)
    U = ift(EE.*ft(U));
    U = gflux(tau,U);
  end
  U = ift(EE.*ft(U));
  U = gflux(tau/2,U);
