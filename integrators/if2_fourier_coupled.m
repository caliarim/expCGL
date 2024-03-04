function U = if2_fourier_coupled(U,M,g,T,m,ft,ift)

  tau = T/m;

  EE{1} = exp(tau*M{1});
  EE{2} = exp(tau*M{2});
  t = 0;
  for jj=1:m
    Gn = g(U);
    Utmp{1} = ift(EE{1}.*ft(U{1}+tau*Gn{1}));
    Utmp{2} = ift(EE{2}.*ft(U{2}+tau*Gn{2}));
    Gn2 = g(Utmp);
    U{1} = ift(EE{1}.*ft(U{1}+(tau/2)*Gn{1}))+(tau/2)*Gn2{1};
    U{2} = ift(EE{2}.*ft(U{2} + (tau/2)*Gn{2}))+(tau/2)*Gn2{2};
    t = t+tau;
  end
