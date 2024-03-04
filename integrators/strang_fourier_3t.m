function U = strang_fourier_3t(U,M,gf1,gf2,T,m,ft,ift)
  tau = T/m;

  EE = exp(tau*M);
  U = gf2(tau/2,U);
  for jj = 1:(m-1)
    U = gf1(tau/2,U);
    U = ift(EE.*ft(U));
    U = gf1(tau/2,U);
    U = gf2(tau,U);
  end
  U = gf1(tau/2,U);
  U = ift(EE.*ft(U));
  U = gf1(tau/2,U);
  U = gf2(tau/2,U);
