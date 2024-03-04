function U = if4_fourier_coupled(U,M,g,T,m,ft,ift)
  tau = T/m;

  EEh{1} = exp((tau/2)*M{1});
  EEh{2} = exp((tau/2)*M{2});

  EE{1} = EEh{1}.*EEh{1};
  EE{2} = EEh{2}.*EEh{2};

  t = 0;
  for jj=1:m
    Gn = g(U);
    Uhat{1} = ft(U{1});
    Uhat{2} = ft(U{2});
    Gnhat{1} = ft(Gn{1});
    Gnhat{2} = ft(Gn{2});
    U2{1} = ift(EEh{1}.*(Uhat{1}+(tau/2)*Gnhat{1}));
    U2{2} = ift(EEh{2}.*(Uhat{2}+(tau/2)*Gnhat{2}));
    Gn2 = g(U2);
    U3{1} = ift(EEh{1}.*Uhat{1})+(tau/2)*Gn2{1};
    U3{2} = ift(EEh{2}.*Uhat{2})+(tau/2)*Gn2{2};
    Gn3 = g(U3);
    U4{1} = ift(EE{1}.*Uhat{1})+tau*ift(EEh{1}.*ft(Gn3{1}));
    U4{2} = ift(EE{2}.*Uhat{2})+tau*ift(EEh{2}.*ft(Gn3{2}));
    Gn4 = g(U4);
    U{1} = ift(EE{1}.*(Uhat{1}+(tau/6)*Gnhat{1}))+(tau/3)*ift(EEh{1}.*ft(Gn2{1}+Gn3{1}))+(tau/6)*Gn4{1};
    U{2} = ift(EE{2}.*(Uhat{2}+(tau/6)*Gnhat{2}))+(tau/3)*ift(EEh{2}.*ft(Gn2{2}+Gn3{2}))+(tau/6)*Gn4{2};
    t = t+tau;
  end
