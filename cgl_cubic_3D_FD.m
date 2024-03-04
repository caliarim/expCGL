clear all
close all

addpath('integrators')
addpath('extern/KronPACK/src/')

d = 3;
n = [128,128,128];
ld = [0,0,0];
rd = [100,100,100];
T = 100;
m = 1000;
method = 'split4'; % strang, if2, split4, if4, rk4

alpha1 = 1;
beta1 = 2;
alpha2 = 1;
alpha3 = -1;
beta3 = 0.2;

for mu = 1:d
  x{mu} = linspace(ld(mu),rd(mu),n(mu)+1).';
  x{mu} = x{mu}(2:n(mu)+1);
  h(mu) = (rd(mu)-ld(mu))/n(mu);
  D2{mu} = spdiags(ones(n(mu),1)*([-1,16,-30,16,-1]/(12*h(mu)^2)),-2:2,n(mu),n(mu));
  D2{mu}(1,1:5) = [-15,-4,14,-6,1]/(12*h(mu)^2);
  D2{mu}(n(mu)-1,(n(mu)-5):n(mu)) = [1,-6,14,-4,-15,10]/(12*h(mu)^2);
  D2{mu}(n(mu),(n(mu)-4):n(mu)) = [1,-8/3,-6,56,-15-100/3]/(12*h(mu)^2);

  M{mu} = full((alpha1+1i*beta1)*D2{mu})+(alpha2/d)*eye(n(mu));
end

[X{1:d}] = ndgrid(x{1:d});

g = @(u) (alpha3+1i*beta3)*(u.*u).*conj(u);
f = @(u) kronsumv(u,M)+g(u);

gflux = @(t,u) exp(-((alpha3+1i*beta3)/(2*alpha3))*log(abs(1-(2*t*alpha3)*(u.*conj(u))))).*u;

rng('default');
U0 = randn(n)/5000;

disp(['Computing solution with ',method,'...'])
tic
switch method
  case 'strang'
    U = strang_FD(U0,M,gflux,T,m);

  case 'if2'
    U = if2_FD(U0,M,g,T,m);

  case 'split4'
    U = split4_FD(U0,M,gflux,T,m);

  case 'if4'
    U = if4_FD(U0,M,g,T,m);

  case 'rk4'
    U = rk4(U0,f,T,m);

  otherwise
    error('Method not known')

end
wc_time = toc;
disp(sprintf('Wall-clock time (seconds): %.2f', wc_time))

figure;
idx = [2,1,3];
isosurface(permute(X{1},idx),permute(X{2},idx),X{3},permute(abs(U),idx),0.3)
xlim([ld(1) rd(1)])
ylim([ld(2) rd(2)])
zlim([ld(3) rd(3)])
xlabel('x_1')
ylabel('x_2')
zlabel('x_3')
title(sprintf('|u| at t=%.2f',T))
drawnow


rmpath('integrators')
rmpath('extern/KronPACK/src/')
