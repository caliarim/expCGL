function U=split4_fourier(U,M,gflux,T,m,ft,ift)
  tau = T/m;
  EEh = exp((tau/2)*M);
  EE = EEh.*EEh;
  for jj = 1:m
    Utmp = gflux(tau/2,ift(EE.*ft(gflux(tau/2,U))));
    U = (4/3)*gflux(tau/4,ift(EEh.*ft(gflux(tau/2,ift(EEh.*ft(gflux(tau/4,U)))))))-...
        (1/3)*Utmp;
  end
