function U=strang_FD(U,M,gflux,T,m)
  d = length(M);
  tau = T/m;
  for mu = 1:d
    E{mu} = expm(tau*M{mu});
  end
  U = gflux(tau/2,U);
  for jj=1:(m-1)
    U = tucker(U,E);
    U = gflux(tau,U);
  end
  U = tucker(U,E);
  U = gflux(tau/2,U);
