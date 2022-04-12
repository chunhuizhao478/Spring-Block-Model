%Read and plot data
clear all; close all; clc

format long
%Import data file
res1 = importdata('res_pbc_0.05_factors.txt' );
res2 = importdata('res_pbc_0.1_factors.txt');
res3 = importdata('res_pbc_0.15_factors.txt' );
% res4 = importdata('res15e6.txt' );

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
%%%%%%res4%%%%%%%%%
% p4     = res4(:,1)/50e6; 
% theta4 = res4(:,2);     
% u4     = res4(:,3);       
% v4     = res4(:,4);    
% psi4   = res4(:,5);  
% t4     = res4(:,6); 
% tau4   = res4(:,7);

%plot figure
%%%%velocity%%%%
%Overview%
f1 = figure(1);
f1.Position = [0 0 2000 1000];
plot(t1,v1,'k-','LineWidth',1.5);hold on
plot(t2,v2,'b-','LineWidth',1.5);hold on
plot(t3,v3,'r-','LineWidth',1.5);hold on
% semilogx(t4,v4,'-','Color','#77AC30','LineWidth',1.5);
title('slider velocity $$\hat{v}$$ v.s. time $$\hat{t}$$ [PBC]','interpreter','latex','FontSize',20)
xlabel('time $$\hat{t}$$','interpreter','latex','FontSize',20)
ylabel('slider velocity $$\hat{v}$$','interpreter','latex','FontSize',20)
xlim([0,2000])
ylim([0, 7000])
xline(0,'m-.',{'Velocity','Perturbation'},FontSize = 16) %Initial Velocity Perturbation
xline(400,'m-.',{'Injection','Start Increase'},FontSize = 16) %Injection starting point
xline(533,'m-.',{'Injection','Reach Maximum'},FontSize = 16) 
xline(1000,'m-.',{'Extraction','Start Decrease'},FontSize = 16) %Injection end point
xline(1133,'m-.',{'Extraction','Reach Minimum'},FontSize = 16)
legend('$\hat{p}_{max} = 0.05$','$\hat{p}_{max} = 0.1$','$\hat{p}_{max} = 0.15$','Interpreter','latex',FontSize = 16)
ax = gca; 
ax.FontSize = 16;

% %Around injection period%
% figure(2);
% plot(t1,v1,'k-','LineWidth',1.5);hold on
% plot(t2,v2,'b-','LineWidth',1.5);hold on
% plot(t3,v3,'r-','LineWidth',1.5);
% title('slider velocity v.s. time [PBC]')
% xlabel('time')
% ylabel('slider velocity')
% xlim([300,2000])
% ylim([0, inf])
% xline(0,'m-.',{'Velocity','Perturbation'}) %Initial Velocity Perturbation
% xline(400,'m-.',{'Injection','Start Point'}) %Injection starting point
% xline(700,'m-.',{'Injection','End Point'}) %Injection end point
% legend('p=1e7 [20%]','p=1.05e7 [21%]','p=1.1e7 [22%]','','','','','','','','','','','')

%Plot intersesmic period%
[pks1,locs1] = findpeaks(v1,'MinPeakDistance',10);
[pks2,locs2] = findpeaks(v2,'MinPeakDistance',10);
[pks3,locs3] = findpeaks(v3,'MinPeakDistance',10);
% [pks4,locs4] = findpeaks(v4,'MinPeakDistance',10);
% plot(t1(locs),pks,'or')
threshold = 0.1;
pk_index1 = find(pks1>threshold); 
pk_index2 = find(pks2>threshold); 
pk_index3 = find(pks3>threshold); 
% pk_index4 = find(pks4>threshold); 

pks1 = pks1(pk_index1);
pks2 = pks2(pk_index2);
pks3 = pks3(pk_index3);
% pks4 = pks4(pk_index4);

locs1 = locs1(pk_index1);
locs2 = locs2(pk_index2);
locs3 = locs3(pk_index3);
% locs4 = locs4(pk_index4);

% plot(t2(locs2),pks2,'or'); 

intersesmic_list1 = diff(t1(locs1));
intersesmic_list2 = diff(t2(locs2));
intersesmic_list3 = diff(t3(locs3));
% intersesmic_list4 = diff(t3(locs4));

t1_plot = t1(locs1); t1_plot = t1_plot(2:end);
t2_plot = t2(locs2); t2_plot = t2_plot(2:end);
t3_plot = t3(locs3); t3_plot = t3_plot(2:end);
% t4_plot = t4(locs4); t4_plot = t4_plot(2:end);

f3 = figure(3);
f3.Position = [0 0 1000 500];
plot1 = loglog(t1_plot,intersesmic_list1,'ks','LineWidth',1.5); hold on
plot2 = loglog(t2_plot,intersesmic_list2,'bo','LineWidth',1.5); hold on
plot3 = loglog(t3_plot,intersesmic_list3,'rd','LineWidth',1.5); hold on
% plot4 = loglog(t4_plot,intersesmic_list4,'^-','Color','#77AC30','LineWidth',1.5);
plot1.MarkerFaceColor = plot1.Color;
plot2.MarkerFaceColor = plot2.Color;
plot3.MarkerFaceColor = plot3.Color;
% plot4.MarkerFaceColor = plot4.Color;
title('interseismic time vs time $$\hat{t}$$ [PBC]','interpreter','latex','FontSize',20)
xlabel('time $$\hat{t}$$','interpreter','latex','FontSize',20)
ylabel('interseismic time','interpreter','latex','FontSize',20)
xlim([400,2000])
ylim([-1, 1400])
xline(0,'m-.',{'Velocity','Perturbation'},FontSize = 16) %Initial Velocity Perturbation
xline(400,'m-.',{'Injection','Start Increase'},FontSize = 16) %Injection starting point
xline(533,'m-.',{'Injection','Reach Maximum'},FontSize = 16) 
xline(1000,'m-.',{'Extraction','Start Decrease'},FontSize = 16) %Injection end point
xline(1133,'m-.',{'Extraction','Reach Minimum'},FontSize = 16)
legend('$\hat{p}_{max} = 0.05$','$\hat{p}_{max} = 0.1$','$\hat{p}_{max} = 0.15$','Interpreter','latex',FontSize = 16)

locs1 = locs1/1e3;
locs2 = locs2/1e3;
locs3 = locs3/1e3;
% locs4 = locs4/1e3 + 400;

f3_2 = figure(32);
f3_2.Position = [0 0 1000 500];
plot1 = semilogy(2:1:length(pks1),intersesmic_list1,'ks-','LineWidth',1.5); hold on
plot2 = semilogy(2:1:length(pks2),intersesmic_list2,'bo-','LineWidth',1.5); hold on
plot3 = semilogy(2:1:length(pks3),intersesmic_list3,'rd-','LineWidth',1.5); hold on
% plot4 = semilogy(1:1:length(pks4),locs4,'^-','Color','#77AC30','LineWidth',1.5);
plot1.MarkerFaceColor = plot1.Color;
plot2.MarkerFaceColor = plot2.Color;
plot3.MarkerFaceColor = plot3.Color;
% plot4.MarkerFaceColor = plot4.Color;
title('event vs interseismic time [PBC]','interpreter','latex','FontSize',20)
xlabel('event number','interpreter','latex','FontSize',20)
ylabel('intersesmic time','interpreter','latex','FontSize',20)
xlim([0, 110])
ylim([0, 600])
legend('$\hat{p}_{max} = 0.05$','$\hat{p}_{max} = 0.1$','$\hat{p}_{max} = 0.15$','Interpreter','latex',FontSize = 16)

f3_3 = figure(39);
f3_3.Position = [0 0 1000 500];
plot1 = semilogy(1:1:length(pks1),locs1,'ks','LineWidth',1.5); hold on
plot2 = semilogy(1:1:length(pks2),locs2,'bo','LineWidth',1.5); hold on
plot3 = semilogy(1:1:length(pks3),locs3,'rd','LineWidth',1.5); hold on
% plot4 = semilogy(1:1:length(pks4),locs4,'^-','Color','#77AC30','LineWidth',1.5);
plot1.MarkerFaceColor = plot1.Color;
plot2.MarkerFaceColor = plot2.Color;
plot3.MarkerFaceColor = plot3.Color;
% plot4.MarkerFaceColor = plot4.Color;
title('event vs time $$\hat{t}$$ [PBC]','interpreter','latex','FontSize',20)
xlabel('event number','interpreter','latex','FontSize',20)
ylabel('time $$\hat{t}$$','interpreter','latex','FontSize',20)
xlim([0, 110])
ylim([0, 2000])
yline(0,'m-.',{'Velocity','Perturbation'},FontSize = 16) %Initial Velocity Perturbation
yline(400,'m-.',{'Injection','Start Increase'},FontSize = 16) %Injection starting point
yline(533,'m-.',{'Injection','Reach Maximum'},FontSize = 16) 
yline(1000,'m-.',{'Extraction','Start Decrease'},FontSize = 16) %Injection end point
yline(1133,'m-.',{'Extraction','Reach Minimum'},FontSize = 16)
legend('$\hat{p}_{max} = 0.05$','$\hat{p}_{max} = 0.1$','$\hat{p}_{max} = 0.15$','Interpreter','latex',FontSize = 16)

% %Plot displacement%
% figure(4);
% plot(t1,u1,'k-','LineWidth',1.5);hold on
% plot(t2,u2,'b-','LineWidth',1.5);hold on
% plot(t3,u3,'r-','LineWidth',1.5);
% title('displacement v.s. time [PBC]')
% xlabel('time')
% ylabel('displacement')
% xlim([0,1000])
% xline(0,'m-.',{'Velocity','Perturbation'}) %Initial Velocity Perturbation
% xline(400,'m-.',{'Injection','Start Point'}) %Injection starting point
% xline(700,'m-.',{'Injection','End Point'}) %Injection end point
% legend('p=1e7 [20%]','p=1.05e7 [21%]','p=1.1e7 [22%]','','','')

% %Plot psi%
% figure(5);
% plot(t1,psi1,'k-','LineWidth',1.5);hold on
% plot(t2,psi2,'b-','LineWidth',1.5);hold on
% plot(t3,psi3,'r-','LineWidth',1.5);hold on
% plot(t4,psi4,'-','Color','#77AC30','LineWidth',1.5);
% title('psi v.s. time [PBC]')
% xlabel('time')
% ylabel('psi')
% xlim([400,2000])
% xline(0,'m-.',{'Velocity','Perturbation'}) %Initial Velocity Perturbation
% xline(400,'m-.',{'Injection','Start Increase'}) %Injection starting point
% xline(533,'m-.',{'Injection','Reach Maximum'}) 
% xline(1000,'m-.',{'Extraction','Start Decrease'}) %Injection end point
% xline(1133,'m-.',{'Extraction','Reach Minimum'})
% yline(0.8,'m--',{'k == k_{cr}'})
% legend('p_{max}=7.5e6 [15%]','p_{max}=1e7 [20%]','p_{max}=1.25e7 [25%]','p_{max}=1.5e7 [30%]','','','')
% ax = gca; 
% ax.FontSize = 16;

% %Plot theta%
% figure(6);
% plot(t1,theta1,'k-','LineWidth',1.5);hold on
% plot(t2,theta2,'b-','LineWidth',1.5);hold on
% plot(t3,theta3,'r-','LineWidth',1.5);
% title('state variable v.s. time [PBC]')
% xlabel('time')
% ylabel('state variable')
% xlim([0,1000])
% xline(0,'m-.',{'Velocity','Perturbation'}) %Initial Velocity Perturbation
% xline(400,'m-.',{'Injection','Start Point'}) %Injection starting point
% xline(700,'m-.',{'Injection','End Point'}) %Injection end point
% legend('p=1e7 [20%]','p=1.05e7 [21%]','p=1.1e7 [22%]','','','')

% %Plot friction force%
% figure(7);
% plot(t1,tau1,'k-','LineWidth',1.5);hold on
% plot(t2,tau2,'b-','LineWidth',1.5);hold on
% plot(t3,tau3,'r-','LineWidth',1.5);
% title('friction force v.s. time [PBC]')
% xlabel('time')
% ylabel('friction force')
% xlim([0,1000])
% xline(0,'m-.',{'Velocity','Perturbation'}) %Initial Velocity Perturbation
% xline(400,'m-.',{'Injection','Start Point'}) %Injection starting point
% xline(700,'m-.',{'Injection','End Point'}) %Injection end point
% legend('p=1e7 [20%]','p=1.05e7 [21%]','p=1.1e7 [22%]','','','')

%Plot pressure%
figure(8);
plot(t1,p1,'k-','LineWidth',1.5);hold on
plot(t2,p2,'b-','LineWidth',1.5);hold on
plot(t3,p3,'r-','LineWidth',1.5);hold on
% plot(t4,p4,'-','Color','#77AC30','LineWidth',1.5);
title('pressure $$\hat{P}$$ v.s. time $\hat{t}$ [PBC]','interpreter','latex','FontSize',20);
xlabel('time $$\hat{t}$$','interpreter','latex','FontSize',20);
ylabel('pressure $$\hat{P}$$','interpreter','latex','FontSize',20);
xlim([0,2000])
ylim([0,0.2])
xline(0,'m-.',{'Velocity','Perturbation'},FontSize = 16) %Initial Velocity Perturbation
xline(400,'m-.',{'Injection','Start Increase'},FontSize = 16) %Injection starting point
xline(533,'m-.',{'Injection','Reach Maximum'},FontSize = 16) 
xline(1000,'m-.',{'Extraction','Start Decrease'},FontSize = 16) %Injection end point
xline(1133,'m-.',{'Extraction','Reach Minimum'},FontSize = 16)
yline(2.5e6/50e6,'m--','5 percent of normal stress, $$k < k_{cr}$$','Interpreter','latex',FontSize = 16)
yline(5e6/50e6,'m--','10 percent of normal stress, $$k = k_{cr}$$','Interpreter','latex',FontSize = 16)
yline(7.5e6/50e6,'m--','15 percent of normal stress, $$k > k_{cr}$$','Interpreter','latex',FontSize = 16)
legend('$\hat{p}_{max} = 0.05$','$\hat{p}_{max} = 0.1$','$\hat{p}_{max} = 0.15$','Interpreter','latex',FontSize = 16)
ax = gca; 
ax.FontSize = 16;

