clear all; clc; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%case1
%%injection duration: t = 400 - 700
%%simulation time to t=1000
%%Introduce memory function d\psi/dt
%main_ode45%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global p_cur %Current pressure
global max_p %Maximum pressure 
%%%%%%%%%%Without Pressure%%%%%%%%%%
%add for loop
max_p_list = [7.5e6]; % 20% 40% 60% of initial normal stress %1.5e7-> 30%
for i = 1:length(max_p_list)
    max_p = max_p_list(i);
    %time
    t_min = 400;
    t_max = 2000;
    t_step = 1e-3; %tuned
    numofval = length(t_min:t_step:t_max)-1;
    %Constant Variables
    V_p = 1e-9; %constant plate velocity
    f_o = 0.6; %friction coefficient
    a = 0.015; b = 0.025; %friction coefficients
    v_o =10^-9; %reference velocity
    L = 10^-6; %characteristic length scale 
    theta_o = L / v_o; %reference state variable
    sigma_initial = 50e6; %normal stress
    p_cur = 0;%initial pressure at center
    % max_p = 5e6; %Maximum pressure during injection period
    p_back = 0;
    sigma_np = sigma_initial - p_cur;
    k_cr=sigma_np*(b-a)/L; %%%k_cr vary with pressure!%%%
    k_s = 0.8*sigma_initial*(b-a)/L; %k_s=0.8*k_cr; %%%k does not vary with pressure%%% system stiffness
    M = k_s;
    %%%%%%Modified dtheta/dt%%%%%%
    alpha = 0.53; %coupling between normal stress and state variable
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Nondimensionlization
    coeff1 = k_s * L ^ 2 / ( M * v_o ^ 2);
    coeff2 = sigma_initial / ( k_s * L );
    coeff3 = V_p / v_o;
    coeff4 = L / ( theta_o * v_o);
    %initial condition
    %v0 = 0 / v_o;
    psi0 = 1;
    v0 = 0 / v_o - 1e-11 / v_o; %initial relative velocity
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
    f2 = @(t,y) [ coeff4 -  abs(y(3)+coeff3) * y(1)
                  y(3);
                  coeff1 * (-y(2) - sign(y(3)+coeff3)* ( sigma_initial / ( k_s * L ) ) * ( f_o + a * log(abs(y(3)+coeff3)) + b * log(abs(y(1)))))];
                  
    %Calculate Jacobian
    %jacobian = @(t,y) [-coeff1*(2*coeff2*dirac(coeff3 + y(1))*(f_o + b*log(abs(y(2))) + a*log(abs(coeff3 + y(1)))) + (a*coeff2*sign(coeff3 + y(1))^2)/abs(coeff3 + y(1))) -(b*coeff1*coeff2*sign(coeff3 + y(1))*sign(y(2)))/abs(y(2)) -coeff1;
    %                         -y(2)*sign(coeff3 + y(1))                 -abs(coeff3 + y(1))         0   ;
    %                         1                           0                 0   ];
    %options = odeset('Jacobian',jacobian);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %define list for store results
    res = zeros(numofval,7);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %pre-simulate to ss (t=400)  %use it when debugging%
    % t_pre_ss = 0;
    % t_new_ss = 400;
    % [t_ss,y_ss] = ode23s(f2,[t_pre_ss,t_new_ss],[theta0,u0,v0]);
    % theta0 = y_ss(end,1); %7.915313563264281
    % u0 = y_ss(end,2); %-62.437688273140900
    % v0 = y_ss(end,3); %-0.999969422154186
    
    % %%%%%%%%%%%%%For k=0.8k_cr%%%%%%%%%%%%%%%%%%%%%%%%%%
    theta0 = 7.915313563264281;
    u0 = -62.437688273140900;
    v0 = -0.999969422154186;
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    psi0 = 1;
    disp('Pre-simulation Done')
    %%%%%%%%%%%%%time integration%%%%%%%%%%%%%%%%%%%%%%%
    p_list = []; t_list=[];
    %initialize time
    t_pre = 400; %Get the time before perform increment
    t_inject_start = 400;
    t_inject_end = 800; %Tanh function, it doesn't mean to stop injection
    %initialize counter
    counter = 1;
    while t_pre < t_max
        %Get pressure at injection point (@t_pre)
        [p_input] = pressure_center(t_pre,t_inject_start,t_inject_end);
        %Get pressure distribution (@t_pre)
        %p0(ceil(end/2)) = p_input;
        %Update diffusion equation (@t_pre -> @t_new=t_pre+t_step)
        %[t1,p] = ode45(f1,[t_pre,t_pre+t_step],p0);
        %Use RHS of diffusion equation to get dp/dt for SBM (@t_pre)
        %f1_rhs = D_c * A_mat_sparse * p0(end,:)';
        %f1_rhs_center = f1_rhs(ceil(end/2),1);
        %Update SBM system of equations (@t_pre -> @t_new=t_pre+t_step)
        p_cur = p_input;
        %y results are for the (@t_new=t_pre+t_step)
        [t2,y] = ode15s(@(t,y)sbm_sys_test1(y,coeff1,coeff2,coeff3,coeff4,sigma_initial,f_o,a,b),[t_pre,t_pre+t_step],[theta0,u0,v0,psi0]);
        %[t2,y] = ode23s(@(t,y)sbm_sys_test1(y,coeff1,coeff3,coeff4,coeff5,sigma_initial,k_s,L,f_o,a,b,f1_rhs_center),[t_pre,t_pre+t_step],[theta0,u0,v0]);
        %Get new pressure at inject location (@t_new=t_pre+t_step)
        %Store Results (@t_new=t_pre+t_step)
        res(counter,1) = p_cur;
        res(counter,2) = y(end,1); %theta
        res(counter,3) = y(end,2); %u
        res(counter,4) = y(end,3); %v
        res(counter,5) = y(end,4); %psi
        res(counter,6) = t_pre+t_step; %time
        res(counter,7) = - sign(y(end,3)+coeff3) * ( coeff2 * y(end,4) ) * ( f_o + a * log(abs(y(end,3)+coeff3)) + b * log(abs(y(end,1))));
        %Plot slider velocity vs time %check
        if mod(counter,1e3*100) == 0 %plot velocity vs time every 500 time steps
            figure(101)
            plot(res(1:counter,6),res(1:counter,4)); hold on
            title('test_pbc'+string(max_p))
            pause(0.5)
            figure(102)
            plot(res(1:counter,6),res(1:counter,1)); hold on
            pause(0.5)
        end
         %Plot pressure spatial distribution
        if mod(counter,1e3*100) == 0 
            p_list = [p_list ; p(end,:)];
            t_list = [t_list t_pre];
            figure(101)
            plot(p(end,:)); hold on
            title('test'+string(i)+'spatial')
            pause(0.5)
        end
%         p0 = p(end,:);
        theta0 = y(end,1);
        u0 = y(end,2);
        v0 = y(end,3);
        psi0 = y(end,4);
        %p_pre = p_new;
        counter = counter + 1;
        t_pre = t_pre + t_step;
        %output time
        disp(t_pre)
    end
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %tau = coeff1 .* sign(res(:,2)+coeff3).* ( (sigma_initial - res(:,1)) ./ ( k_s .* L ) ) .* ( f_o + a .* log(abs(res(:,2)+coeff3)) + b .* log(abs(res(:,3))));
    
    %%%%%%%%%%%%%%post-processing%%%%%%%%%%%%%%%%%%%%
    %save res
    fileID_1 = fopen('res'+string(max_p)+'test_ct'+'.txt','w');
    [rows1, columns1] = size(res);
    for row = 1 : rows1
        for col = 1 : columns1
	        fprintf(fileID_1, '%f ', res(row, col));
        end
        fprintf(fileID_1, '\n');
    end
    fclose(fileID_1);
end


print('All Done!')