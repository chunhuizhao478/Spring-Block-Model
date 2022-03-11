%-----------------------------%
%Start from physical quantites
%-----------------------------%
%Solve 1D diffusion equation
% clear all; clc;
function [dx,err] = diffusion_main(n_num)
%Generate 1D mesh
n_start = -200; %m
n_end = 200; %m
n_num = n_num; %m
x = linspace(n_start,n_end,n_num); %m
dx = x(2) - x(1);
%Constant variable (for now)
phi = 0.01; %porosity 
beta = 6.4e-10; %fluid compressibility %[ct:=sum of both solid and fluid compressibility] %Pa^-1 1Pa = 1kg.m^-1.s^-2
mu = 1e-3; %Pa s -> Pa h %beta * mu -> Pa * Pa^-1, no need to scale
D_c = 1 / (phi * beta * mu); %coefficient s^-1 -> h^-1
p0 = zeros(n_num-2,1); 
source = zeros(size(p0)); 
m_input = 3e-6;
rho = 1e3;
%Q_input = 6e-2; 
%Obtain A matrix (pure diffusion)
%-------------------------------%
[A_mat] = Amat(x,D_c);
A_mat_sparse = sparse(A_mat);
%-------------------------------%
%Obtain A matrix (during injection, include flow rate)
%---------------------------------------------------%
[A_mat_injection] = Amat_injection(x,D_c);
A_mat_sparse_injection = sparse(A_mat_injection);
%---------------------------------------------------%
p_list = []; t_list=[]; err_list = [];
%initialize time
t_min = 0; %s 
t_max = 1; %s 
t_step = 0.001; %s CFL = dt / dx^2 < 1 dx=1 dt = 0.1
numofval = length(t_min:t_step:t_max)-1;
t_pre = 0; %s  %Get the time before perform increment
t_inject_start = 0; %s
t_inject_end = 50; %s 
res = zeros(numofval,2);
counter = 1;
%load analytical solution
k_mid =  1.28e-14;
t_compare = 1;
D = k_mid / (phi * beta * mu);
S = m_input / ( rho * phi * beta );
x_real = linspace(0,n_end,(n_num-1)/2+1); x_real = x_real(1:end-1);
while t_pre < t_max
        %update pressure at the center
        if t_pre > t_inject_start && t_pre < t_inject_end
            source(ceil(end/2)) = m_input / ( rho * phi * beta ) * 1/dx ;
            A_mat_current = A_mat_sparse_injection;
        else
            source(ceil(end/2)) = 0;
            A_mat_current = A_mat_sparse;
        end
        options = odeset('RelTol',1e-12,'AbsTol',1e-12);
        f1 = @(t,p) A_mat_current * p + source;
        %Update diffusion equation (@t_pre -> @t_new=t_pre+t_step)
        [t1,p] = ode45(f1,[t_pre,t_pre+t_step],p0,options);
        %get pressure at the injection point
        p_cur = p(end,ceil(end/2));
        %store the pressure spatial distribution
        t_now = t_pre + t_step;
        if counter == round(t_compare/t_step)
%              %Flip the data
%              p_flip = flip(p(end,:));
%              %Assemble together
%              p_spatialtotal = horzcat(p_flip(1:end-1),p(end,:));
             p_spatialtotal = p(end,:);
             %Store the results
             p_list = [p_list ; p(end,:)];
             t_list = [t_list ; t_now];
             %Analytical Solution
             real_p = test_analytical(D,t_compare,x_real,S);
             %Compute Error
             err = (sum(((p_spatialtotal-real_p)/real_p).^2)*dx)^0.5;
             disp('err:')
             disp(err)
             disp('dx:')
             disp(dx)
             %Plot the results
             figure(50)
             real_xcoord = linspace(-200,200,n_num-2);
             fd_xcoord = linspace(-200,200,n_num-2);
             plot(fd_xcoord, p_spatialtotal,'k-'); hold on
             plot(real_xcoord, real_p,'r-');
             title('Diffusion Equation w Source Term')
             xlabel('x coordinate')
             ylabel('pressure')
             legend('Finite Difference','Analytical Solution')
        end
        %store the pressure at the injection point
         res(counter,1) = t_pre+t_step;
         res(counter,2) = p_cur;
%          if mod(counter,10) == 0
%              figure(101)
%              plot(res(1:counter,1),res(1:counter,2)); hold on
%              pause(0.5)
%          end
%          figure(9)
%          plot(err_list)
         %Update parameters
         p0 = p(end,:);
         counter = counter + 1;
         t_pre = t_pre + t_step;
         %output time
%          disp(t_pre)
end
% figure()
% plot(res(1:counter-1,1),res(1:counter-1,2));
% title('pressure at center')
% xlabel('time')
% ylabel('pressure')
writematrix(p_spatialtotal,'p_spatial_out.txt','WriteMode','append')
writematrix(res,'res_out.txt','WriteMode','append')
end