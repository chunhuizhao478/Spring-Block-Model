function [f2] = sbm_sys_test1(y,coeff1,coeff2,coeff3,coeff4,sigma_initial,f_o,a,b)
%unknown [theta,u,v,psi]
global p_cur 
f2 = [ coeff4 -  abs(y(3)+coeff3) * y(1); 
       y(3);
       coeff1 * (-y(2) - sign(y(3)+coeff3) * ( coeff2 * y(4) ) * ( f_o + a * log(abs(y(3)+coeff3)) + b * log(abs(y(1)))))
       -0.1 * abs(y(3)+coeff3) * ( y(4) - (sigma_initial - p_cur) / (sigma_initial) )];
end

%coeff1 * (-y(2) - sign(y(3)+coeff3)* ( coeff2 * y(4) ) * ( f_o + a * log(abs(y(3)+coeff3)) + b * log(abs(y(1)))))