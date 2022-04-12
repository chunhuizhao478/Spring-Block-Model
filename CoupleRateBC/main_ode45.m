clear all; clc; close all;
format long
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%main_ode45%
%Diffusion equation w source term
%Include theta/pressure coupling
dt_factor = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global p_cur p_pre %Current pressure
% global max_p %Maximum pressure 
%%%%%%%%%%Without Pressure%%%%%%%%%%
%add for loop
injection_end_time_list = [];
m_input_list = [0.75 0.5]; %pmax=1.5e7 2.4e7 3.3e7  (sigma_initial=5.0e7)
% Q_input_list = [1e4];
for i = 1:length(m_input_list)
    m_input = m_input_list(i);
    %time
    t_min = 1000;
    t_max = 1e4;
    t_step = 1e-3; %tunedp %1e-3
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
    p_pre = 0;
    % max_p = 5e6; %Maximum pressure during injection period
    p_back = 0;
    sigma_np = sigma_initial - sigma_initial * p_cur;
    k_cr=sigma_np*(b-a)/L; %%%k_cr vary with pressure!%%%
    k_s = 0.9*sigma_initial*(b-a)/L; %k_s=0.9*k_cr; %%%k does not vary with pressure%%% system stiffness
    M = k_s;
    rho = 1e3; %kg/m^3
    m_o = 3e-6; %kg/s 3e-7 %res3em6
    %%%%%%Modified dtheta/dt%%%%%%
    alpha = 0.53; %coupling between normal stress and state variable
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Nondimensionlization
    coeff1 = k_s * L ^ 2 / ( M * v_o ^ 2);
    coeff2 = 1 / ( k_s * L );
    coeff3 = V_p / v_o;
    coeff4 = L / ( theta_o * v_o);
    %initial condition
    %v0 = 0 / v_o;
    psi0 = 1;
    v0 = 0 / v_o - 1e-11 / v_o; %initial relative velocity
    theta0 = coeff4 / abs(v0 + coeff3); %initial state variable
    F_tau0 = sigma_initial * coeff2 * ( f_o + a * log(abs(v0+coeff3)) + b * log(abs(theta0)));%initial fricition force
    u0 = -F_tau0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%Include Pressure%%%%%%%%%%
    %Generate 1D mesh
    %------------------%
    %Increase domain%
    n_start = -200 * 2;
    n_end = 200 * 2;
    n_num = 400 * 2 + 1;
    %------------------%
    x = linspace(n_start,n_end,n_num);
    dx = x(2) - x(1);
    %Constant variable (for now)
    phi = 0.01; %porosity
    beta = 6.4e-10; %fluid compressibility %[ct:=sum of both solid and fluid compressibility]
    mu = 1e-3; %Pa s
    D_c = 1 / (phi * beta * mu); %coefficient
    %Modified coefficient of diffusivity
    Lx = 100 / 2; %characteristic length scale, take it as half of the damage zone size
    D_c_modifier = L / v_o;
    D_c = D_c_modifier * D_c;
    %Obtain A matrix
%     [A_mat] = Amat(x,D_c);
%     A_mat_sparse = sparse(A_mat);
    %source term
    p0 = zeros(n_num-2,1); 
    source = zeros(size(p0)); 
    %Q = 10/3 * 10^4; %Assumed %%%test cases: 2/3 * 10^4(2e6) || 4/3 * 10^4(4e6) || 6/3 * 10^4(6e6) || 8/3 * 10^4(8e6) || 10/3 * 10^4(10e6)   
    %source(ceil(end/2)) = Q_input;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%Combine Together%%%%%%%%%%
    %Define system of equations
    %%diffusion equation%%
    %     f1 = @(t,p) D_c * A_mat_sparse * p;
    
    %%SBM system of equations%% 
    %%unknown: [u,theta,v] %%
    %no pressure/theta coupling%
%     f2 = @(t,y) [ coeff4 -  abs(y(3)+coeff3) * y(1) + ( alpha * sigma_initial * y(1) ) / ( b * (sigma_initial - p_cur * sigma_initial) ) * ( p_cur - p_pre) / t_step; 
%                   y(3);
%                   coeff1 * (-y(2) - sign(y(3)+coeff3) * ( coeff2 * (sigma_initial - p_cur * sigma_initial) ) * ( f_o + a * log(abs(y(3)+coeff3)) + b * log(abs(y(1))))); ];
                  
    %Calculate Jacobian
    %jacobian = @(t,y) [-coeff1*(2*coeff2*dirac(coeff3 + y(1))*(f_o + b*log(abs(y(2))) + a*log(abs(coeff3 + y(1)))) + (a*coeff2*sign(coeff3 + y(1))^2)/abs(coeff3 + y(1))) -(b*coeff1*coeff2*sign(coeff3 + y(1))*sign(y(2)))/abs(y(2)) -coeff1;
    %                         -y(2)*sign(coeff3 + y(1))                 -abs(coeff3 + y(1))         0   ;
    %                         1                           0                 0   ];
    %options = odeset('Jacobian',jacobian);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %define list for store results
    res = zeros(round(numofval*dt_factor),7);
    
    %%%%%%%%%%%%Build A matrix before loop%%%%%%%%%%%%%%
    %injection
    %Nondimensionalized factor
    [A_mat_injection] = Amat_injection(x,D_c,Lx);
    A_mat_sparse_injection = sparse(A_mat_injection);
    %before injection
%     [A_mat] = Amat(x,D_c,Lx);
%     A_mat_sparse = sparse(A_mat);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %pre-simulate to ss (t=400)  %use it when debugging%
%     t_pre_ss = 0;
%     t_new_ss = 400;
%     [t_ss,y_ss] = ode23s(f2,[t_pre_ss,t_new_ss],[theta0,u0,v0]);
%     theta0 = y_ss(end,1); %7.915313563264281
%     u0 = y_ss(end,2); %-62.437688273140900
%     v0 = y_ss(end,3); %-0.999969422154186
    
    % %%%%%%%%%%%%%For k=0.9k_cr @ t=1000 %%%%%%%%%%%%%%%%%%%%%%%%%%
      theta0 = 6.075725949937613;
      u0 = -57.383006716665875;
      v0 = -0.999811679466282;
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%     psi0 = 1;
    disp('Pre-simulation Done')
    %%%%%%%%%%%%%time integration%%%%%%%%%%%%%%%%%%%%%%%
    p_list = []; t_list=[]; %store spatial distribution
    %initialize time
    t_pre = 1000; %Get the time before perform increment
    t_inject_start = 1000;
    %t_inject_end = 3000;
    %initialize counter
    counter = 1;
    %Set threshold value % 7.5e6 -> 15% normal stress % 1e7 -> 20% normal stress
    threshold_value = 7.5e6;
    injection_tag = true;
    while t_pre < t_max
        %update pressure at the center
        if t_pre > t_inject_start && p_cur * sigma_initial <= threshold_value && injection_tag == true
             source(ceil(end/2)) = L * m_o / (sigma_initial * v_o * phi * beta * rho) * m_input * 1/dx;
             A_mat_current = A_mat_sparse_injection;
        else
             source(ceil(end/2)) = 0;
             A_mat_current = A_mat_sparse_injection;
        end
        if p_cur * sigma_initial > threshold_value
           injection_tag = false; %Can't access again
           injection_end_time = t_pre;
           injection_end_time_list = [injection_end_time_list t_pre];
        end
        f1 = @(t,p) A_mat_current * p + source;
        %Update diffusion equation (@t_pre -> @t_new=t_pre+t_step)
        [t1,p] = ode45(f1,[t_pre,t_pre+t_step],p0);
        %Update SBM system of equations (@t_pre -> @t_new=t_pre+t_step)
        p_cur = p(ceil(end/2));
        %y results are for the (@t_new=t_pre+t_step)
        %options = odeset('RelTol',1e-6,'AbsTol',1e-6);
        %[t2,y] = ode15s(@(t,y)sbm_sys_test1(y,coeff1,coeff2,coeff3,coeff4,sigma_initial,f_o,a,b),[t_pre,t_pre+t_step],[theta0,u0,v0,psi0],options);
        [t2,y] = ode15s(@(t,y)sbm_sys_test1(y,coeff1,coeff2,coeff3,coeff4,sigma_initial,f_o,a,b,t_step,alpha),[t_pre,t_pre+t_step],[theta0,u0,v0]);
        %Get new pressure at inject location (@t_new=t_pre+t_step)
        %Store Results (@t_new=t_pre+t_step)
        t_new = t_pre + t_step;

        res(counter,1) = p_cur;
        res(counter,2) = y(end,1); %theta
        res(counter,3) = y(end,2); %u
        res(counter,4) = y(end,3); %v
%         res(counter,5) = y(end,4); %psi
        res(counter,6) = t_pre+t_step; %time
%         res(counter,7) = - sign(y(end,3)+coeff3) * ( coeff2 * y(end,4) ) * ( f_o + a * log(abs(y(end,3)+coeff3)) + b * log(abs(y(end,1))));

        %Plot slider velocity vs time %check
        if mod(counter,round(1/t_step)*200) == 0 %plot velocity vs time every 500 time steps
            figure(99 + 199 * i)
            plot(res(1:counter,6),res(1:counter,4)); hold on
            title('test'+string(i)+'velocity')
            pause(0.5)
            disp('#-------------------#')
            disp('# Current Velocity  #')
            disp('#-------------------#')
            disp(y(end,3))
        end
        if mod(counter,round(1/t_step)*200) == 0 %plot pressure vs time every 500 time steps
            figure(100 + 199 * i)
            plot(res(1:counter,6),res(1:counter,1)); hold on
            title('test'+string(i)+'pressure_center')
            pause(0.5)
        end
%         if mod(counter,1e3*20) == 0 %plot psi vs time every 500 time steps
%             figure(101 + 199 * i)
%             plot(res(1:counter,6),res(1:counter,5)); hold on
%             title('test'+string(i)+'psi')
%             pause(0.5)
%         end
        %Plot pressure spatial distribution
        if mod(counter,round(1/t_step)*200) == 0 
            p_list = [p_list ; p(end,:)];
            t_list = [t_list t_pre+t_step];
            figure(101)
            plot(p(end,:)); hold on
            title('test'+string(i)+'spatial')
            pause(0.5)
        end
        p0 = p(end,:);
        theta0 = y(end,1);
        u0 = y(end,2);
        v0 = y(end,3);
%         psi0 = y(end,4);
        counter = counter + 1;
        t_pre = t_pre + t_step;
        p_pre = p_cur;
        
        %output time
        disp(t_pre)
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%post-processing%%%%%%%%%%%%%%%%%%%
    %save results
    %save results %must be double quote1
    writematrix(t_list,"spatial"+"_"+"fixpmax"+"_"+string(m_input)+"_"+"time.txt")
    writematrix(p_list,"spatial"+"_"+"fixpmax"+"_"+string(m_input)+"_"+"pressure.txt")
    writematrix(res,"res"+"_"+"fixpmax"+"_"+string(m_input)+"_"+"factors.txt")
end

%print injection end time
disp('####Injection_End_Time_List####')
disp(injection_end_time_list)

print('All Done!')

