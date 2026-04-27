clear all
close all
clc;
load worktemp_l.mat
delta = 1-(1-0.05)^(1/4);
alpha = 1/3; 
ypc = worktemp_l.mled(:,2);
xpc = worktemp_l.mled(:,3);
hpc = worktemp_l.mled(:,4);


% annualizing the variables

for i=4:size(ypc,1)
    ypc_a(i) = sum(ypc(i-3:i))';
    xpc_a(i) = sum(xpc(i-3:i))';
end

ypc_a=ypc_a';
xpc_a=xpc_a';


kss = mean(xpc)/delta;
kpc = ypc*NaN;
kpc(1) = kss;
for i = 1:size(ypc,1)-1
    kpc(i+1) = (1-delta)*kpc(i)+xpc(i);
end

dypc_a = diff(log(ypc_a));
dkpc = diff(log(kpc));
dhpc = diff(log(hpc));

dA = dypc_a-alpha*dkpc-(1-alpha)*dhpc;
% plot(xpc)
% plot(ypc)
% plot(hpc)

t = 1980.25:0.25:2015;
t0 = find(t==2008.25);

plot(t(t0:end)',dypc_a(t0-1:end)); hold on;
plot(t(t0:end)',dA(t0-1:end)); hold on;
plot(t(t0:end)',dkpc(t0-1:end)); hold on;
plot(t(t0:end)',dhpc(t0-1:end)); hold off;
legend('\Delta Y','\Delta A','\Delta K','\Delta H');
title('Growth accounting for the US, 2000-2015');


%%% o chari quer os dados sem ser per capita e sem ser detrended pela tech 