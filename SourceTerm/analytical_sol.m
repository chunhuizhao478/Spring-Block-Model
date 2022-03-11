function [p_real] = analytical_sol(D,t,x,S)
%Analytical solution of diffusion equation (homogeneous permeability)
p_real = S * ( sqrt( t / ( D * pi ) ) * exp( - x.^2 / ( 4 * D * t ) ) - 1/(2*D).* x .* erfc( x./sqrt(4*D*t)) );
end