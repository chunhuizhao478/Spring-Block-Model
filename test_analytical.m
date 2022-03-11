function [p_spatialtotal] = test_analytical(D,Q,k,mu,t,x)

p_real = analytical_sol(D,Q,k,mu,t,x);
%Flip the data
p_flip = flip(p_real);
%Assemble together
p_spatialtotal = horzcat(p_flip(1:end-1),p_real);
