function U=rk2_coupled(U,f,T,m,alpha)
  tau = T/m;
  t = 0;
  for jj=1:m
    k1 = f(U);
    Utmp{1} = U{1}+(alpha*tau)*k1{1};
    Utmp{2} = U{2}+(alpha*tau)*k1{2};
    k2 = f(Utmp);
    U{1} = U{1}+tau*((1-1/(2*alpha))*k1{1}+(1/(2*alpha))*k2{1});
    U{2} = U{2}+tau*((1-1/(2*alpha))*k1{2}+(1/(2*alpha))*k2{2});
    t = t+tau;
  end
