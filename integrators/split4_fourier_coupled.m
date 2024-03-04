function U = split4_fourier_coupled(U,M,g,T,m,ft,ift)
  tau = T/m;

  EEh{1} = exp(tau/2*M{1});
  EEh{2} = exp(tau/2*M{2});
  
  EE{1} = EEh{1}.*EEh{1};
  EE{2} = EEh{2}.*EEh{2};
  t = 0;
  for jj=1:m
    Utmp = rk4_coupled(U,g,tau/2,1);
    Utmp{1} = ift(EE{1}.*ft(Utmp{1}));
    Utmp{2} = ift(EE{2}.*ft(Utmp{2}));
    Utmp = rk4_coupled(Utmp,g,tau/2,1);

    U = rk4_coupled(U,g,tau/4,1);
    U{1} = ift(EEh{1}.*ft(U{1}));
    U{2} = ift(EEh{2}.*ft(U{2}));
    U = rk4_coupled(U,g,tau/2,1);
    U{1} = ift(EEh{1}.*ft(U{1}));
    U{2} = ift(EEh{2}.*ft(U{2}));
    U = rk4_coupled(U,g,tau/4,1);
    U{1} = (4/3)*U{1}-(1/3)*Utmp{1};
    U{2} = (4/3)*U{2}-(1/3)*Utmp{2};
    t = t+tau;
  end
