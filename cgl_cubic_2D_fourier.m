clear all
close all

addpath('integrators')

d = 2;

n = [256,256];
ld = [0,0];
rd = [100,100];
T = 30;
m = 50;
method = 'split4'; % strang, if2, split4, if4, rk2, rk4

alpha1 = 1;
beta1 = 2;
alpha2 = 1;
alpha3 = -1;
beta3 = 0.2;

for mu = 1:d
  x{mu} = linspace(ld(mu),rd(mu),n(mu)+1).';
  x{mu} = x{mu}(1:n(mu));
  h(mu) = (rd(mu)-ld(mu))/(n(mu));
  lambda{mu} = 1i*2*pi*(-n(mu)/2:n(mu)/2-1)'/(rd(mu)-ld(mu));
  lambda2{mu} = lambda{mu}.*lambda{mu};
end

[X{1:d}] = ndgrid(x{1:d});
[L2{1:d}] = ndgrid(lambda2{1:d});

M = (alpha1+1i*beta1)*(L2{1}+L2{2})+alpha2;

nc = sqrt(prod(rd-ld))/prod(n);
ft = @(u) fftshift(fft2(u))*nc;
ift = @(u) ifft2((ifftshift(u)))/nc;

g = @(u) (alpha3+1i*beta3)*((u.*u).*conj(u));
f = @(u) ift(M.*ft(u))+g(u);

gflux = @(t,u) exp(-((alpha3+1i*beta3)/(2*alpha3))*log(abs(1-(2*t*alpha3)*(u.*conj(u))))).*u;

rng('default');
U0 = randn(n)/5000;

disp(['Computing solution with ',method,'...'])
tic
switch method
  case 'strang'
    U = strang_fourier(U0,M,gflux,T,m,ft,ift);

  case 'if2'
    U = if2_fourier(U0,M,g,T,m,ft,ift);

  case 'split4'
    U = split4_fourier(U0,M,gflux,T,m,ft,ift);

  case 'if4'
    U = if4_fourier(U0,M,g,T,m,ft,ift);

  case 'rk2'
    U = rk2(U0,f,T,m,1);

  case 'rk4'
    U = rk4(U0,f,T,m);

  otherwise
    error('Method not known')

end
wc_time = toc;
disp(sprintf('Wall-clock time (seconds): %.2f', wc_time))

figure;
pcolor(X{1},X{2},abs(U))
shading interp
colorbar
title(sprintf('|u| at t=%.2f',T))
xlabel('x_1')
ylabel('x_2')
drawnow

rmpath('integrators')
