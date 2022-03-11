%-----------------------------%
%Start from physical quantites
%-----------------------------%
%Solve 1D diffusion equation
clear all; clc;
%Generate 1D mesh
n_start = -200; %m
n_end = 200; %m
n_num = 401; %m
x = linspace(n_start,n_end,n_num); %m
%Constant variable (for now)
phi = 0.01; %porosity 
beta = 6.4e-10; %fluid compressibility %[ct:=sum of both solid and fluid compressibility] %Pa^-1 1Pa = 1kg.m^-1.s^-2
mu = 1e-3 / 3600; %Pa s -> Pa h %beta * mu -> Pa * Pa^-1, no need to scale
D_c = 1 / (phi * beta * mu); %coefficient s^-1 -> h^-1
p0 = zeros(n_num-2,1); 
source = zeros(size(p0)); 
Q_input = 6e-2 * 3600; %m/s -> m/h 20bpm - 1bpm = 0.16m^3/min
%Obtain A matrix (pure diffusion)
%-------------------------------%
[A_mat] = Amat(x,D_c);
A_mat_sparse = sparse(A_mat);
%-------------------------------%
%Obtain A matrix (during injection, include flow rate)
%---------------------------------------------------%
[A_mat_injection] = Amat_injection(x,mu,D_c);
A_mat_sparse_injection = sparse(A_mat_injection);
%---------------------------------------------------%
p_list = []; t_list=[];
%initialize time
t_min = 0; %s -> h  400 * L/v_o = 4e5
t_max = 1e6; %s -> h
t_step = 1; %s CFL = dt / dx^2 < 1 dx=1 dt = 0.1
numofval = length(t_min:t_step:t_max)-1;
t_pre = 0; %s -> h %Get the time before perform increment
t_inject_start = 0; %s -> h
t_inject_end = 5.6e4; %s -> h
res = zeros(numofval,2);
counter = 1;
while t_pre < t_max
        %update pressure at the center
        if t_pre > t_inject_start && t_pre < t_inject_end
            source(ceil(end/2)) = Q_input;
            A_mat_current = A_mat_sparse_injection;
        else
            source(ceil(end/2)) = 0;
            A_mat_current = A_mat_sparse;
        end
        f1 = @(t,p) A_mat_current * p + source;
        %Update diffusion equation (@t_pre -> @t_new=t_pre+t_step)
        [t1,p] = ode45(f1,[t_pre,t_pre+t_step],p0);
        %get pressure at the injection point
        p_cur = p(end,ceil(end/2));
        %store the pressure spatial distribution
        if mod(counter,1) == 0 
             p_list = [p_list ; p(end,:)];
%              t_list = [t_list t_pre];
%              figure(100)
%              plot(p(end,:)); hold on
%              pause(0.5)
        end
        %store the pressure at the injection point
        res(counter,1) = t_pre+t_step;
        res(counter,2) = p_cur;
        if mod(counter,1e3) == 0
            figure(101)
            plot(res(1:counter,1),res(1:counter,2)); hold on
            pause(0.5)
        end
        %Update parameters
        p0 = p(end,:);
        counter = counter + 1;
        t_pre = t_pre + t_step;
        %output time
        disp(t_pre)
end
%save res
fileID_1 = fopen('p_input'+string(Q_input)+'.txt','w');
[rows1, columns1] = size(res);
for row = 1 : rows1
    for col = 1 : columns1
        fprintf(fileID_1, '%f ', res(row, col));
    end
    fprintf(fileID_1, '\n');
end
fclose(fileID_1);
writematrix(p_list,'pressure_spatial'+string(Q_input)+'.txt')