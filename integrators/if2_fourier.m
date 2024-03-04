function U=if2_fourier(U,M,g,T,m,ft,ift)
  tau = T/m;
  t = 0;
  EE = exp(tau*M);
  for jj=1:m
    Gn = g(U);
    Utmp = ift(EE.*ft(U+tau*Gn));
    U = ift(EE.*ft(U+(tau/2)*Gn))+(tau/2)*g(Utmp);
    t = t+tau;
  end
