function U=split4_fourier_3t(U,M,gf1,gf2,T,m,ft,ift)
  tau = T/m;
  EEh = exp((tau/2)*M);
  EE = EEh.*EEh;
  for jj = 1:m
    Utmp=gf2(tau/2,gf1(tau/2,ift(EE.*ft(gf1(tau/2,gf2(tau/2,U))))));
    U=(4/3)*gf2(tau/4,gf1(tau/4,ift(EEh.*ft(gf1(tau/4,gf2(tau/2,gf1(tau/4,ift(EEh.*ft(gf1(tau/4,gf2(tau/4,U)))))))))))-...
    (1/3)*Utmp;
  end
