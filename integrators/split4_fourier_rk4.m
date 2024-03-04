function U=split4_fourier_rk4(U,M,g,T,m,ft,ift)
  tau = T/m;
  EEh = exp((tau/2)*M);
  EE = EEh.*EEh;
  t = 0;
  for jj = 1:m
    Utmp=rk4(ift(EE.*ft(rk4(U,g,tau/2,1))),g,tau/2,1);
    U= (4/3)*rk4(ift(EEh.*ft(rk4(ift(EEh.*ft(rk4(U,g,tau/4,1))),g,tau/2,1))),g,tau/4,1)-...
    (1/3)*Utmp;
    t = t + tau;
  end
