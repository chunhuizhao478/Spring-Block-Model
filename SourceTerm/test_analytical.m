function [p_spatialtotal] = test_analytical(D,t,x,S)

% %Generate 1D mesh
% n_start = 0; %m
% n_end = 200; %m
% n_num = 201; %m
% x = linspace(n_start,n_end,n_num); %m

p_real = analytical_sol(D,t,x,S);
%Flip the data
p_flip = flip(p_real);
%Assemble together
p_spatialtotal = horzcat(p_flip(1:end-1),p_real);

% plot(p_spatialtotal)
