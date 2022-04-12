%Read and plot data
clear all; close all; clc

format long
%Import data file
res1 = importdata("res_fixpmax_1_factors.txt" );
res2 = importdata("res_fixpmax_0.75_factors.txt");
res3 = importdata("res_fixpmax_0.5_factors.txt" );

%Label data set
%%%%%%res1%%%%%%%%%
p1     = res1(:,1); 
theta1 = res1(:,2);     
u1     = res1(:,3);       
v1     = res1(:,4);    
psi1   = res1(:,5);  
t1     = res1(:,6); 
tau1   = res1(:,7);
%%%%%%res2%%%%%%%%%
p2     = res2(:,1); 
theta2 = res2(:,2);     
u2     = res2(:,3);       
v2     = res2(:,4);    
psi2   = res2(:,5);  
t2     = res2(:,6); 
tau2   = res2(:,7);
%%%%%%res3%%%%%%%%%
p3     = res3(:,1); 
theta3 = res3(:,2);     
u3     = res3(:,3);       
v3     = res3(:,4);    
psi3   = res3(:,5);  
t3     = res3(:,6); 
tau3   = res3(:,7);

%plot figure
% 5e4 8e4 1e5
% k b r
%%%%velocity%%%%
%Overview%
f1 = figure(1);
f1.Position = [0 0 1000 500];
plot(t1,v1,'k-','LineWidth',0.5);hold on
plot(t2,v2,'b-','LineWidth',0.5);hold on
plot(t3,v3,'r-','LineWidth',0.5);
title('slider velocity $$\hat{\bar{v}}$$ v.s. time $$\hat{t}$$ [RateBC]','interpreter','latex','FontSize',20)
xlabel('time $$\hat{t}$$','interpreter','latex','FontSize',20)
ylabel('slider velocity $$\hat{\bar{v}}$$','interpreter','latex','FontSize',20)
xlim([0,1e4])
ylim([0, 7000])
xline(1000,'m-.') %Injection starting point
xline(1643.56,'m-.',{'$\hat{m} = 1.0$','End Point'},'Interpreter','latex',FontSize = 16) %Injection end point
xline(2143.97,'m-.',{'$\hat{m} = 0.75$','End Point'},'Interpreter','latex',FontSize = 16) %Injection end point
xline(3567.72,'m-.',{'$\hat{m} = 0.5$','End Point'},'Interpreter','latex',FontSize = 16) %Injection end point
legend('$\hat{m} = 1.0$','$\hat{m} = 0.75$','$\hat{m} = 0.5$','Interpreter','latex',FontSize = 16)
ax = gca; 
ax.FontSize = 16;

% %Around injection period%

%Plot intersesmic period%
[pks1,locs1] = findpeaks(v1,'MinPeakDistance',10);
[pks2,locs2] = findpeaks(v2,'MinPeakDistance',10);
[pks3,locs3] = findpeaks(v3,'MinPeakDistance',10);
% plot(t1(locs),pks,'or')
threshold = 1.0;
pk_index1 = find(pks1>threshold); 
pk_index2 = find(pks2>threshold); 
pk_index3 = find(pks3>threshold); 

pks1 = pks1(pk_index1);
pks2 = pks2(pk_index2);
pks3 = pks3(pk_index3);

locs1 = locs1(pk_index1);
locs2 = locs2(pk_index2);
locs3 = locs3(pk_index3);

% plot(t2(locs2),pks2,'or'); 

intersesmic_list1 = diff(t1(locs1));
intersesmic_list2 = diff(t2(locs2));
intersesmic_list3 = diff(t3(locs3));

t1_plot = t1(locs1); t1_plot = t1_plot(2:end);
t2_plot = t2(locs2); t2_plot = t2_plot(2:end);
t3_plot = t3(locs3); t3_plot = t3_plot(2:end);

f3=figure(3);
f3.Position = [0 0 1000 500];
loglog(t1_plot,intersesmic_list1,'ko-','LineWidth',1.5); hold on
loglog(t2_plot,intersesmic_list2,'bo-','LineWidth',1.5); hold on
loglog(t3_plot,intersesmic_list3,'ro-','LineWidth',1.5); 
title('interseismic time vs time $$\hat{t}$$ [RateBC]','interpreter','latex','FontSize',20)
xlabel('time $$\hat{t}$$','interpreter','latex','FontSize',20)
ylabel('interseismic time','interpreter','latex','FontSize',20)
xlim([1000,7000])
%ylim([0, inf])
xline(1000,'m-.') %Injection starting point
xline(1643.56,'m-.',{'$\hat{m} = 1.0$','End Point'},'Interpreter','latex',FontSize = 16) %Injection end point
xline(2143.97,'m-.',{'$\hat{m} = 0.75$','End Point'},'Interpreter','latex',FontSize = 16) %Injection end point
xline(3567.72,'m-.',{'$\hat{m} = 0.5$','End Point'},'Interpreter','latex',FontSize = 16) %Injection end point
legend('$\hat{m} = 1.0$','$\hat{m} = 0.75$','$\hat{m} = 0.5$','Interpreter','latex',FontSize = 16)
ax = gca; 
ax.FontSize = 16;

locs1 = locs1/1e3 ;
locs2 = locs2/1e3 ;
locs3 = locs3/1e3 ;

f3_2 = figure(32);
f3_2.Position = [0 0 1000 500];
plot1 = semilogy(2:1:length(pks1),intersesmic_list1,'k.-','LineWidth',1.5); hold on
plot2 = semilogy(2:1:length(pks2),intersesmic_list2,'b.-','LineWidth',1.5); hold on
plot3 = semilogy(2:1:length(pks3),intersesmic_list3,'r.-','LineWidth',1.5); 
plot1.MarkerFaceColor = plot1.Color;
plot2.MarkerFaceColor = plot2.Color;
plot3.MarkerFaceColor = plot3.Color;
title('event vs interseismic time [RateBC]','interpreter','latex','FontSize',20)
xlabel('event number','interpreter','latex','FontSize',20)
ylabel('intersesmic time','interpreter','latex','FontSize',20)
xlim([0,500])
ylim([0,510])
% yline(1000,'m-.',{'Injection','Start Point'}) %Injection starting point
% yline(1643.56,'m-.',{'Injection m=1','End Point'}) %Injection end point
% yline(2143.97,'m-.',{'Injection m=0.75','End Point'}) %Injection end point
% yline(3567.72,'m-.',{'Injection m=0.5','End Point'}) %Injection end point
legend('$\hat{m} = 1.0$','$\hat{m} = 0.75$','$\hat{m} = 0.5$','Interpreter','latex',FontSize = 16)
ax = gca; 
ax.FontSize = 16;

f3_3 = figure(33);
f3_3.Position = [0 0 1000 500];
plot1 = semilogy(1:1:length(pks1),locs1,'k.-','LineWidth',1.5); hold on
plot2 = semilogy(1:1:length(pks2),locs2,'b.-','LineWidth',1.5); hold on
plot3 = semilogy(1:1:length(pks3),locs3,'r.-','LineWidth',1.5); 
plot1.MarkerFaceColor = plot1.Color;
plot2.MarkerFaceColor = plot2.Color;
plot3.MarkerFaceColor = plot3.Color;
title('event vs time $$\hat{t}$$ [RateBC]','interpreter','latex','FontSize',20)
xlabel('event number','interpreter','latex','FontSize',20)
ylabel('time $$\hat{t}$$','interpreter','latex','FontSize',20)
xlim([0,inf])
ylim([0,inf])
yline(1000,'m-.',{'Injection','Start Point'}) %Injection starting point
yline(1643.56,'m-.',{'$\hat{m} = 1.0$','End Point'},'Interpreter','latex',FontSize = 16) %Injection end point
yline(2143.97,'m-.',{'$\hat{m} = 0.75$','End Point'},'Interpreter','latex',FontSize = 16) %Injection end point
yline(3567.72,'m-.',{'$\hat{m} = 0.5$','End Point'},'Interpreter','latex',FontSize = 16) %Injection end point
legend('$\hat{m} = 1.0$','$\hat{m} = 0.75$','$\hat{m} = 0.5$','Interpreter','latex',FontSize = 16)
ax = gca; 
ax.FontSize = 16;

% %Plot displacement%
% figure(4);
% loglog(t1,u1,'k-','LineWidth',1.5);hold on
% loglog(t2,u2,'b-','LineWidth',1.5);hold on
% loglog(t3,u3,'r-','LineWidth',1.5);
% title('displacement v.s. time [RateBC]')
% xlabel('time')
% ylabel('displacement')
% xlim([0,2500])
% xline(400,'m-.',{'Injection','Start Point'}) %Injection starting point
% xline(2187,'m-.',{'Injection m=1','End Point'}) %Injection end point
% xline(3555,'m-.',{'Injection m=0.75','End Point'}) %Injection end point
% xline(6930,'m-.',{'Injection m=0.5','End Point'}) %Injection end point
% legend('Q=5e4','Q=8e4','Q=1e5','','','')

% %Plot psi%
% figure(5);
% plot(t1,psi1,'k-','LineWidth',1.5);hold on
% plot(t2,psi2,'b-','LineWidth',1.5);hold on
% plot(t3,psi3,'r-','LineWidth',1.5);
% title('psi v.s. time [RateBC]')
% xlabel('time')
% ylabel('psi')
% xlim([0,1.5e4])
% xline(400,'m-.',{'Injection','Start Point'}) %Injection starting point
% xline(2187,'m-.',{'Injection m=1','End Point'}) %Injection end point
% xline(3555,'m-.',{'Injection m=0.75','End Point'}) %Injection end point
% xline(6930,'m-.',{'Injection m=0.5','End Point'}) %Injection end point
% yline(0.8,'m--',{'k == k_{cr}'})
% legend('$$\hat{m}$$=1.0','$$\hat{m}$$=0.75','$$\hat{m}$$=0.5','Interpreter','latex')
% ax = gca; 
% ax.FontSize = 16;

% %Plot theta%
% figure(6);
% plot(t1,theta1,'k-','LineWidth',1.5);hold on
% plot(t2,theta2,'b-','LineWidth',1.5);hold on
% plot(t3,theta3,'r-','LineWidth',1.5);
% title('state variable v.s. time [RateBC]')
% xlabel('time')
% ylabel('state variable')
% xlim([0,2500])
% xline(0,'m-.',{'Velocity','Perturbation'}) %Initial Velocity Perturbation
% xline(400,'m-.',{'Injection','Start Point'}) %Injection starting point
% xline(700,'m-.',{'Injection','End Point'}) %Injection end point
% legend('Q=5e4','Q=8e4','Q=1e5','','','')

% %Plot friction force%
% figure(7);
% plot(t1,tau1,'k-','LineWidth',1.5);hold on
% plot(t2,tau2,'b-','LineWidth',1.5);hold on
% plot(t3,tau3,'r-','LineWidth',1.5);
% title('friction force v.s. time [RateBC]')
% xlabel('time')
% ylabel('friction force')
% xlim([0,2500])
% xline(0,'m-.',{'Velocity','Perturbation'}) %Initial Velocity Perturbation
% xline(400,'m-.',{'Injection','Start Point'}) %Injection starting point
% xline(700,'m-.',{'Injection','End Point'}) %Injection end point
% legend('Q=5e4','Q=8e4','Q=1e5','','','')

%Plot pressure%
figure(8);
semilogx(t1,p1,'k-','LineWidth',1.5);hold on
semilogx(t2,p2,'b-','LineWidth',1.5);hold on
semilogx(t3,p3,'r-','LineWidth',1.5);
title('pressure $$\hat{P}$$ v.s. time $\hat{t}$ [RateBC]','interpreter','latex','FontSize',20);
xlabel('time $$\hat{t}$$','interpreter','latex','FontSize',20);
ylabel('pressure $$\hat{P}$$','interpreter','latex','FontSize',20);
xlim([1000,1e4])
ylim([0,0.2])
%xline(1000,'m-.',{'Injection','Start Point'},'Interpreter','latex',FontSize = 16) %Injection starting point
xline(1643.56,'m-.',{'$\hat{m} = 1.0$','End Point'},'Interpreter','latex',FontSize = 16) %Injection end point
xline(2143.97,'m-.',{'$\hat{m} = 0.75$','End Point'},'Interpreter','latex',FontSize = 16) %Injection end point
xline(3567.72,'m-.',{'$\hat{m} = 0.5$','End Point'},'Interpreter','latex',FontSize = 16) %Injection end point
yline(0.05,'m--','5 percent of normal stress, $$k = k_{cr}$$','Interpreter','latex',FontSize = 16)
yline(0.1,'m--','10 percent of normal stress, $$k = k_{cr}$$','Interpreter','latex',FontSize = 16)
yline(0.15,'m--','15 percent of normal stress, $$k > k_{cr}$$','Interpreter','latex',FontSize = 16)
legend('$\hat{m} = 1.0$','$\hat{m} = 0.75$','$\hat{m} = 0.5$','Interpreter','latex',FontSize = 16)
ax = gca; 
ax.FontSize = 16;

