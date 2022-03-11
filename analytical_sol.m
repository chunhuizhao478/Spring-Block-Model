function [p_real] = analytical_sol(D,Q,k,mu,t,x)
%Analytical solution of diffusion equation (homogeneous permeability)
p_real = (Q * mu / k) * ( sqrt( 4 * D * t / pi ) * exp( -x.^2 / ( 4 * D * t ) ) - x .* erfc( x ./ sqrt( 4 * D * t ) ) );
end