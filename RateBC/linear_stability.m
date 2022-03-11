%Linear Stability
syms L theta_o v_o v V_p theta k M u f_o a b c psi sigma_np sigma_n
theta_dot = L / ( theta_o * v_o ) - abs( v + V_p / v_o ) * theta;
u_dot = v;
v_dot = ( k * L ^ 2 ) / ( M * v_o ^ 2 ) * ( - u - ( psi * sigma_n ) / ( k * L ) * ( f_o + a * log( v + V_p / v_o ) + b * log( theta )));
psi_dot = - 1 / c * ( v + V_p / v_o ) * ( psi - sigma_np / sigma_n );

f_x_prime = [ diff(theta_dot, theta) diff(theta_dot, u) diff(theta_dot, v) diff(theta_dot, psi); 
              diff(u_dot, theta)     diff(u_dot, u)     diff(u_dot,v)      diff(u_dot, psi); 
              diff(v_dot, theta)     diff(v_dot, u)     diff(v_dot,v)      diff(v_dot, psi); 
              diff(psi_dot, theta)   diff(psi_dot, u)   diff(psi_dot,v)    diff(psi_dot, psi)]

%Example
syms p_dot lambda t Coeff_V Coeff_Theta alpha p kappa
eqn1 = v_o / a * (b + (alpha * p_dot) / ( v_o * (sigma_n - p)) - kappa / ( sigma_n - p)) * Coeff_V + v_o ^ 2 / a * Coeff_Theta - lambda * Coeff_V;
eqn2 = - ( b + (alpha * p_dot) / ( v_o * (sigma_n - p))) * Coeff_V - v_o * Coeff_Theta - lambda * Coeff_Theta;
% eqn1 = subs(eqn1,Coeff_V,1);
% eqn1 = subs(eqn1,Coeff_Theta,1);
% eqn2 = subs(eqn2,Coeff_V,1);
% eqn2 = subs(eqn2,Coeff_Theta,1);
eqn = simplify(eqn1 + eqn2);

(Coeff_V*v_o*(b + kappa/(p - sigma_n) - (alpha*p_dot)/(v_o*(p - sigma_n))))/a - Coeff_V*lambda - (Coeff_V*v_o^2*(b - (alpha*p_dot)/(v_o*(p - sigma_n))))/(a*(lambda + v_o));