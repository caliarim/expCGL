function U=if4_fourier(U,M,g,T,m,ft,ift)
  tau = T/m;
  EEh = exp((tau/2)*M);
  EE = EEh.*EEh;
  t = 0;
  for jj=1:m
    Uhat = ft(U);
    Gn = g(U);
    Gnhat = ft(Gn);
    U2 = ift(EEh.*(Uhat+(tau/2)*Gnhat));
    Gn2 = g(U2);
    U3 = ift(EEh.*Uhat)+(tau/2)*Gn2;
    Gn3 = g(U3);
    U4 = ift(EE.*Uhat)+tau*ift(EEh.*ft(Gn3));
    U = ift(EE.*(Uhat+(tau/6)*Gnhat))+(tau/3)*ift(EEh.*ft(Gn2+Gn3))+(tau/6)*g(U4);
    t = t+tau;
  end
