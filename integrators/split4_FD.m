function U=split4_FD(U,M,gflux,T,m)
  d = length(M);
  tau = T/m;
  t = 0;
  for mu = 1:d
    Eh{mu} = expm((tau/2)*M{mu});
    E{mu} = Eh{mu}*Eh{mu};
  end
  for jj=1:m
    Utmp = gflux(tau/2,tucker(gflux(tau/2,U),E));
    U = (4/3)*gflux(tau/4,tucker(gflux(tau/2,tucker(gflux(tau/4,U),Eh)),Eh))-(1/3)*Utmp;
    t = t+tau;
  end
