function U = strang_fourier_coupled(U,M,g,T,m,ft,ift)
  tau = T/m;

  EE{1} = exp(tau*M{1});
  EE{2} = exp(tau*M{2});
  U = rk4_coupled(U,g,tau/2,1);
  for jj = 1:(m-1)
    U{1} = ift(EE{1}.*ft(U{1}));
    U{2} = ift(EE{2}.*ft(U{2}));
    U = rk4_coupled(U,g,tau,1);
  end
  U{1} = ift(EE{1}.*ft(U{1}));
  U{2} = ift(EE{2}.*ft(U{2}));
  U = rk4_coupled(U,g,tau/2,1);
