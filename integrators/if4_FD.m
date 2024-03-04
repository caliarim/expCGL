function U=if4_FD(U,M,g,T,m)
  d = length(M);
  tau = T/m;
  t = 0;
  for mu = 1:d
    Eh{mu} = expm((tau/2)*M{mu});
    E{mu} = Eh{mu}*Eh{mu};
  end
  for jj=1:m
    Gn = g(U);
    U2 = tucker(U+(tau/2)*Gn,Eh);
    Gn2 = g(U2);
    U3 = tucker(U,Eh)+(tau/2)*Gn2;
    Gn3 = g(U3);
    U4 = tucker(U,E)+tau*tucker(Gn3,Eh);
    U = tucker(U+(tau/6)*Gn,E)+(tau/3)*tucker(Gn2+Gn3,Eh)+(tau/6)*g(U4);
    t = t + tau;
  end
