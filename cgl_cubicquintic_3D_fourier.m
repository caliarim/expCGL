clear all
close all

addpath('integrators')

d = 3;

n = [128,128,128];
ld = [-12,-12,-12];
rd = [12,12,12];
T = 20;
m = 80;
method = 'if4'; % strang, strang_3t, if2, split4, split4_3t, if4, rk4

alpha1 = 1/2;
beta1 = 1/2;
alpha2 = -1/2;
alpha3 = 2.52;
beta3 = 1;
alpha4 = -1;
beta4 = -0.11;

for mu = 1:d
  x{mu} = linspace(ld(mu),rd(mu),n(mu)+1).';
  x{mu} = x{mu}(1:n(mu));
  h(mu) = (rd(mu)-ld(mu))/(n(mu));
  lambda{mu} = 1i*2*pi*(-n(mu)/2:n(mu)/2-1)'/(rd(mu)-ld(mu));
  lambda2{mu} = lambda{mu}.*lambda{mu};
end

[X{1:d}] = ndgrid(x{1:d});
[L2{1:d}] = ndgrid(lambda2{1:d});

M = (alpha1+1i*beta1)*(L2{1}+L2{2}+L2{3})+alpha2;

nc = sqrt(prod(rd-ld))/prod(n);
ft = @(u) fftshift(fftn(u))*nc;
ift = @(u) ifftn(ifftshift(u))/nc;

g = @(u) gfun_cubicquintic(u,alpha3,beta3,alpha4,beta4);
f = @(u) ift(M.*ft(u))+g(u);

gflux1 = @(t,u) exp(-((alpha3+1i*beta3)/(2*alpha3))*log(abs(1-(2*t*alpha3)*(u.*conj(u))))).*u;
gflux2 = @(t,u) exp(-((alpha4+1i*beta4)/(4*alpha4))*log(abs(1-(4*t*alpha4)*((u.*conj(u)).*(u.*conj(u)))))).*u;

delta = 1.2;
rho0 = 6;
omega = 2.5;
eta = 5;
kappa = 3;

rho = sqrt(X{1}.^2 + X{2}.^2);
theta = atan2(X{2},X{1});
U0 = delta*sech(sqrt((rho-rho0).^2+X{3}.^2)/omega).*cos(eta*theta).*exp(1i*kappa*theta);

t = 0;
figure;
idx = [2,1,3];
isosurface(permute(X{1},idx),permute(X{2},idx),X{3},permute(abs(U0),idx),sqrt(0.5))
xlim([ld(1) rd(1)])
ylim([ld(2) rd(2)])
zlim([ld(3) rd(3)])
xlabel('x_1')
ylabel('x_2')
zlabel('x_3')
title('|u| at t=0')
drawnow

disp(['Computing solution with ',method,'...'])
tic
switch method
  case 'strang'
    U = strang_fourier_rk4(U0,M,g,T,m,ft,ift);

  case 'strang_3t'
    U = strang_fourier_3t(U0,M,gflux1,gflux2,T,m,ft,ift);

  case 'if2'
    U = if2_fourier(U0,M,g,T,m,ft,ift);

  case 'split4'
    U = split4_fourier_rk4(U0,M,g,T,m,ft,ift);

  case 'split4_3t'
    U = split4_fourier_3t(U0,M,gflux1,gflux2,T,m,ft,ift);

  case 'if4'
    U = if4_fourier(U0,M,g,T,m,ft,ift);

  case 'rk4'
    U = rk4(U0,f,T,m);

  otherwise
    error('Method not known')

end
wc_time = toc;
disp(sprintf('Wall-clock time (seconds): %.2f', wc_time))

figure;
idx = [2,1,3];
isosurface(permute(X{1},idx),permute(X{2},idx),X{3},permute(abs(U),idx),sqrt(0.5))
xlim([ld(1) rd(1)])
ylim([ld(2) rd(2)])
zlim([ld(3) rd(3)])
xlabel('x_1')
ylabel('x_2')
zlabel('x_3')
title(sprintf('|u| at t=%.2f',T))
drawnow

rmpath('integrators')
