function U=rk4_coupled(U,f,T,m)
  tau = T/m;
  t = 0;
  for jj=1:m
    k1 = f(U);
    Utmp{1} = U{1}+(tau/2)*k1{1};
    Utmp{2} = U{2}+(tau/2)*k1{2};
    k2 = f(Utmp);
    Utmp{1} = U{1}+(tau/2)*k2{1};
    Utmp{2} = U{2}+(tau/2)*k2{2};
    k3 = f(Utmp);
    Utmp{1} = U{1}+tau*k3{1};
    Utmp{2} = U{2}+tau*k3{2};
    k4 = f(Utmp);
    U{1} = U{1} + (tau/6)*(k1{1}+2*k2{1}+2*k3{1}+k4{1});
    U{2} = U{2} + (tau/6)*(k1{2}+2*k2{2}+2*k3{2}+k4{2});
    t = t + tau;
  end
