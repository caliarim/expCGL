function U = strang_fourier_rk4(U,M,g,T,m,ft,ift)
  tau = T/m;

  EE = exp(tau*M);
  U = rk4(U,g,tau/2,1);
  for jj = 1:(m-1)
    U = ift(EE.*ft(U));
    U = rk4(U,g,tau,1);
  end
  U = ift(EE.*ft(U));
  U = rk4(U,g,tau/2,1);
