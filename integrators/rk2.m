function U=rk2(U,f,T,m,alpha)
  tau = T/m;
  t = 0;
  for jj=1:m
    k1 = f(U);
    k2 = f(U+(alpha*tau)*k1);
    U = U+tau*((1-1/(2*alpha))*k1+(1/(2*alpha))*k2);
    t = t+tau;
  end
