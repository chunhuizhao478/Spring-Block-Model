clear all; clc; %close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%case1
%%injection duration: t = 400 - 700
%%simulation time to t=1000
%%Introduce memory function d\psi/dt
%main_ode45%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global p_cur
%%%%%%%%%%Without Pressure%%%%%%%%%%
%time
t_min = 0;
t_max = 400;
t_step = 1e-3; %tuned
numofval = length(t_min:t_step:t_max)-1;
%Constant Variables
V_p = 1e-9; %constant plate velocity
f_o = 0.6; %friction coefficient
a = 0.015; b = 0.025; %friction coefficients
v_o =1e-9; %reference velocity
L = 1e-6; %characteristic length scale 
theta_o = L / v_o; %reference state variable
sigma_initial = 51e6; %normal stress
p_cur = 0;%initial pressure at center
p_back = 0;
sigma_np = sigma_initial - p_cur;
k_cr=sigma_np*(b-a)/L; %%%k_cr is also varying with pressure!%%%
k_s=0.8*k_cr;
M =1*k_s;
%%%%%%Modified dtheta/dt%%%%%%
alpha = 0.53; %coupling between normal stress and state variable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Nondimensionlization
coeff1 = k_s * L ^ 2 / ( M * v_o ^ 2);
coeff2 = ( sigma_initial - p_cur ) / ( k_s * L );
coeff3 = V_p / v_o;
coeff4 = L / ( theta_o * v_o);
coeff5 = alpha * L / ( b * v_o );
%initial condition
%v0 = 0 / v_o;
psi0 = 1;
v0 = 0 / v_o - 1e-10 / v_o; %initial relative velocity
theta0 = coeff4 / abs(v0 + coeff3); %initial state variable
F_tau0 = coeff2 * ( f_o + a * log(abs(v0+coeff3)) + b * log(abs(theta0)));%initial fricition force
u0 = -F_tau0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%Include Pressure%%%%%%%%%%
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
p0 = zeros(n_num-2,1); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%Combine Together%%%%%%%%%%
%Define system of equations
%%diffusion equation%%
f1 = @(t,p) D_c * A_mat_sparse * p;

%%SBM system of equations%% 
%%unknown: [u,theta,v] %%
%no pressure/theta coupling%
f2 = @(t,y)   [ 
              
              coeff4 -  abs(y(3)+coeff3) * y(1); 
              y(3);
              coeff1 * (-y(2) - sign(y(3)+coeff3)* ( sigma_initial / ( k_s * L ) ) * ( f_o + a * log(abs(y(3)+coeff3)) + b * log(abs(y(1)))))
             
              ];
              
              %-0.1 * abs(y(3)+coeff3) * ( y(4) - (sigma_initial - p_cur) / (sigma_initial) )];

res = zeros(numofval,5);

%options = odeset('RelTol',1e-1,'AbsTol',1e-2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%pre-simulate to ss (t=400)  %use it when debugging%
t_pre_ss = 0;
t_new_ss = 400; 

dt = 1e-6;
tspan = t_pre_ss : dt : t_new_ss;

%tspan = [t_pre_ss t_new_ss];

[t_ss,y_ss] = ode23s(f2,tspan,[theta0,u0,v0]);
theta0 = y_ss(end,1);
u0 = y_ss(end,2);
v0 = y_ss(end,3);
psi0 = 1;
disp('Pre-simulation Done')