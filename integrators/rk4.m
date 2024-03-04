function U=rk4(U,f,T,m)
  tau = T/m;
  t = 0;
  for jj=1:m
    k1 = f(U);
    k2 = f(U+(tau/2)*k1);
    k3 = f(U+(tau/2)*k2);
    k4 = f(U+tau*k3);
    U = U+(tau/6)*(k1+2*k2+2*k3+k4);
    t = t+tau;
  end
