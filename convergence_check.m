clear all; clc; close all;
% n_input = linspace(200,200*30,30) + 1;
n_input = [201];
dx_list = [];
err_list = [];
for item = 1:length(n_input)
    [dx,err] = diffusion_main(n_input(item));
    dx_list = [dx_list;dx];
    err_list = [err_list;err];
    disp('this is step:')
    disp(item)
end
%Plot number of nodes vs err in log-log plot
figure();
n_input_plot = linspace(200,200*30,30) .* 2 + 1;
loglog(n_input_plot,err_list,'-o')
title('Convergence Check')
xlabel('Num of nodes')
ylabel('Relative Error')