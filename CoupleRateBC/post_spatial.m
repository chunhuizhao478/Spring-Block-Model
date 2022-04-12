%Read and plot data
clear all; close all; clc

format long
%Import data file
time = importdata('spatial_case1_time.txt');
pre1 = importdata('spatial_case1_pressure.txt' );
pre2 = importdata('spatial_case2_pressure.txt' );
pre3 = importdata('spatial_case3_pressure.txt' );

x_coord = linspace(-200,200,401);
x_coord = x_coord(2:end-1);

figure();
plot(x_coord,pre1(9,:),'k-','LineWidth',0.5);hold on
plot(x_coord,pre2(11,:),'b-','LineWidth',0.5);hold on
plot(x_coord,pre3(18,:),'r-','LineWidth',0.5);
title('Spatial Pressure Distribution')
xlabel('x coordinate')
ylabel('Pressure')
legend('m=1 @ t=1800','m=0.75 @ t=2200','m=0.5 @ t=3600')
ax = gca; 
ax.FontSize = 16;
