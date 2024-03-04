function U=if2_FD(U,M,g,T,m)
  d = length(M);
  tau = T/m;
  t = 0;
  for mu = 1:d
    E{mu} = expm(tau*M{mu});
  end
  for jj=1:m
    Gn = g(U);
    Utmp = tucker(U+tau*Gn,E);
    U = tucker(U+(tau/2)*Gn,E)+(tau/2)*g(Utmp);
    t = t+tau;
  end
