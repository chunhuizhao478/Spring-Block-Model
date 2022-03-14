%Read and plot data
format long
res = importdata('res.txt');

time = res(:,1);
p_center = res(:,2);
theta = res(:,3);
u = res(:,4);
v = res(:,5);
psi = res(:,6);
tau = res(:,7);

figure(3)
subplot(5,1,1)
plot(time,v)
title('slider velocity v.s. time')
subplot(5,1,2)
plot(time,theta)
title('state variable v.s. time')
subplot(5,1,3)
plot(time,u)
title('relative displacement v.s. time')
subplot(5,1,4)
plot(time,psi)
title('psi v.s. time')
subplot(5,1,5)
plot(time,tau)
title('friction force v.s. time')

figure(4)
plot(time,p_center)
title('injected pressure magitude v.s. time')