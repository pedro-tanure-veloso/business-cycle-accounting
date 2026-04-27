clear all
close all
clc;
load worktemp.mat
delta = 1-(1-0.05)^(1/4);
alpha = 1/3; 
ypc = exp(worktemp.Y(:,1));
xpc = exp(worktemp.Y(:,2));
hpc = exp(worktemp.Y(:,3));


kss = mean(xpc)/delta;
kpc = ypc*NaN;
kpc(1) = kss;
for i = 1:size(ypc,1)-1
    kpc(i+1) = (1-delta)*kpc(i)+xpc(i);
end

dypc = diff(log(ypc));
dkpc = diff(log(kpc));
dhpc = diff(log(hpc));

dA = dypc-alpha*dkpc-(1-alpha)*dhpc;
% plot(xpc)
% plot(ypc)
% plot(hpc)

t = 1980.25:0.25:2015;
t0 = find(t==2008.25);

plot(t(t0:end)',dypc(t0-1:end)); hold on;
plot(t(t0:end)',dA(t0-1:end)); hold on;
plot(t(t0:end)',dkpc(t0-1:end)); hold on;
plot(t(t0:end)',dhpc(t0-1:end)); hold off;
legend('\Delta Y','\Delta A','\Delta K','\Delta H');
title('Growth accounting for the US, 2000-2015');


%%% o chari quer os dados sem ser per capita e sem ser detrended pela tech 