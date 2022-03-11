%Solve 1D diffusion equation
clear all; clc;
%Generate 1D mesh
n_start = -200;
n_end = 200;
n_num = 401;
x = linspace(n_start,n_end,n_num);
%Obtain A matrix
[A_mat] = Amat(x);
A_mat_sparse = sparse(A_mat);
%Constant variable (for now)
phi = 0.01; %porosity
beta = 6.4e-10; %fluid compressibility %[ct:=sum of both solid and fluid compressibility]
mu = 1e-3; %Pa s
D_c = 1 / (phi * beta * mu); %coefficient
%source term
p0 = zeros(n_num-2,1); p0(ceil(end/2)) = 2e6;
%use ode45
%%%define equation
t_min = 0;
t_max = 10;
f = @(t,p) D_c * A_mat_sparse * p;
[t,p] = ode45(f,[t_min,t_max],p0);

figure(1);
plot(x(2:end-1),p(5e2,:),'b-')
hold on
plot(x(2:end-1),p(1e3,:),'m-')
hold on
plot(x(2:end-1),p(1e4,:),'r-')
hold on
plot(x(2:end-1),p(1e5,:),'g-')
hold on
plot(x(2:end-1),p(1e6,:),'k-')
legend('t=2e-4','t=5e-4','t=6e-3','t=6e-2','t=6e-1')
grid on
xlabel('x')
ylabel('pressure')
title('pressure distribution w.r.t time')