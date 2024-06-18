clear all
close all

addpath('integrators')

d = 2;

n = [700,350];
ld = [0,0];
rd = [70,35];
T_st = 25;
T = 15;
m = 1000;
method = 'if4'; % strang, if2, split4, if4, rk2, rk4

alpha0 = -0.4;
alpha1 = 0.125;
beta1 = 0.5;
alpha2 = -0.9;
alpha3 = 1;
beta3 = 0.8;
alpha4 = -0.1;
beta4 = -0.6;
alpha5 = 0.5;

for mu = 1:d
  x{mu} = linspace(ld(mu),rd(mu),n(mu)+1).';
  x{mu} = x{mu}(1:n(mu));
  h(mu) = (rd(mu)-ld(mu))/n(mu);
  lambda{mu} = 1i*2*pi*(-n(mu)/2:n(mu)/2-1)'/(rd(mu)-ld(mu));
  lambda2{mu} = lambda{mu}.*lambda{mu};
end

[X{1:d}] = ndgrid(x{1:d});
[L1{1:d}] = ndgrid(lambda{1:d});
[L2{1:d}] = ndgrid(lambda2{1:d});

% Compute asymptotic solitons in 1D along x direction
M_st = alpha2+(alpha1+1i*beta1)*lambda2{1};

nc_1 = sqrt(rd(1)-ld(1))/n(1);
ft_1 = @(u) fftshift(fft(u))*nc_1;
ift_1 = @(u) ifft(ifftshift(u))/nc_1;

delta = 2.25;
chi = (rd(1)/2-rd(1)/4);
sigma = 2.5;

u0_st = delta*exp(-((x{1}-chi).^2)/(2*(sigma*sigma)));

gfun_st = @(u) gfun_cubicquintic(u,alpha3,beta3,alpha4,beta4);

m_st = 10000;
disp('Computing asymptotic 1D states...')
tic
u_st = if4_fourier(u0_st,M_st,gfun_st,T_st,m_st,ft_1,ift_1);
wc_time_st=toc;
disp('Asymptotic 1D states computed!')
disp(sprintf('Wall-clock time (seconds) asymptotic 1D state: %.2f', wc_time_st))

figure;
plot(x{1},abs(u0_st),'-ro')
hold on
plot(x{1},abs(u_st),'-bx')
xlabel('x_1')
ylabel('|w|')
legend('t=0', sprintf('t=%.2f',T_st))
drawnow

% Now simulate the coupling of quasi-1D states

U0{1} = repmat(u_st,1,n(2));
U0{2} = repmat(flip(u_st),1,n(2));

M{1} = alpha0*L1{1}+(alpha1+1i*beta1)*(L2{1}+L2{2})+alpha2;
M{2} = -alpha0*L1{1}+(alpha1+1i*beta1)*(L2{1}+L2{2})+alpha2;

nc = sqrt(prod(rd-ld))/prod(n);
ft = @(u) fftshift(fft2(u))*nc;
ift = @(u) ifft2(ifftshift(u))/nc;

gfun = @(u) gfun_coupled(u,alpha3,beta3,alpha4,beta4,alpha5);
ffun = @(u) ffun_coupled(u,M,alpha3,beta3,alpha4,beta4,alpha5,ft,ift);

figure
subplot(1,2,1)
pcolor(X{1},X{2},abs(U0{1}))
shading interp
colorbar
xlabel('x_1')
ylabel('x_2')
title('|u| at t=0')
subplot(1,2,2)
pcolor(X{1},X{2},abs(U0{2}))
shading interp
colorbar
xlabel('x_1')
ylabel('x_2')
title('|v| at t=0')
drawnow

disp(['Computing solution with ',method,'...'])
tic
switch method
  case 'strang'
    U = strang_fourier_coupled(U0,M,gfun,T,m,ft,ift);

  case 'if2'
    U = if2_fourier_coupled(U0,M,gfun,T,m,ft,ift);

  case 'split4'
    U = split4_fourier_coupled(U0,M,gfun,T,m,ft,ift);

  case 'if4'
    U = if4_fourier_coupled(U0,M,gfun,T,m,ft,ift);

  case 'rk2'
    U = rk2_coupled(U0,ffun,T,m,1);

  case 'rk4'
    U = rk4_coupled(U0,ffun,T,m);

  otherwise
    error('Method not known')

end
wc_time = toc;
disp(sprintf('Wall-clock time (seconds) simulation: %.2f', wc_time))

figure
subplot(1,2,1)
pcolor(X{1},X{2},abs(U{1}))
shading interp
colorbar
xlabel('x_1')
ylabel('x_2')
title(sprintf('|u| at t=%.2f',T))
subplot(1,2,2)
pcolor(X{1},X{2},abs(U{2}))
shading interp
colorbar
xlabel('x_1')
ylabel('x_2')
title(sprintf('|v| at t=%.2f',T))
drawnow

rmpath('integrators')
