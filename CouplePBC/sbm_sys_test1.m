function [f2] = sbm_sys_test1(y,coeff1,coeff2,coeff3,coeff4,sigma_initial,f_o,a,b,t_step,alpha,Dp_input)
%unknown [theta,u,v]
global p_cur p_pre
f2 = [ coeff4 -  abs(y(3)+coeff3) * y(1) + ( alpha * sigma_initial * y(1) ) / ( b * (sigma_initial - p_cur * sigma_initial) ) * Dp_input; 
       y(3);
       coeff1 * (-y(2) - sign(y(3)+coeff3) * ( coeff2 * (sigma_initial - p_cur * sigma_initial) ) * ( f_o + a * log(abs(y(3)+coeff3)) + b * log(abs(y(1))))) ];
       %-0.1 * abs(y(3)+coeff3) * ( y(4) - (sigma_initial - p_cur) / (sigma_initial) )];
end

%coeff1 * (-y(2) - sign(y(3)+coeff3)* ( coeff2 * y(4) ) * ( f_o + a * log(abs(y(3)+coeff3)) + b * log(abs(y(1)))))

%p_cur is nondimensionized quantity