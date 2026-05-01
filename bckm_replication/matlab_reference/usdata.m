%USDATA    Set up data file for maximum likelihood estimation,
%          after adjusting NIPA output to exclude sales tax and 
%          include consumer durable services.  Store the matrix 
%          `mled' in uszvarq.dat.

%          Ellen McGrattan, 3-15-04
%          Revised, ERM, 12-28-14

clear all;
close all;
clc;

load data/nipa115.dat
load data/nipa116.dat
load data/nipa119.dat
load data/nipa32.dat
load data/nipa33.dat
load data/nipa394.dat
load data/nipa395.dat
load data/atab10d.dat
load data/btab100d.dat
load data/hours.dat
load data/civpop.dat
load data/armed.dat

T     = 1:182;  % 1969:1-2014:2
rGDP  = nipa116(3,T)';
rPCE  = nipa116(4,T)';
pCD   = nipa119(6,T)';
rCD   = nipa115(6,T)'./nipa119(6,T)'*100;
rCND  = nipa115(7,T)'./nipa119(7,T)'*100;
rCS   = nipa115(8,T)'./nipa119(8,T)'*100;
rGPDI = nipa116(9,T)';
rEX   = nipa116(18,T)';
rIM   = nipa116(21,T)';
rG    = nipa116(24,T)';
rGC   = nipa395(4,T)'./nipa394(4,T)'*100;
rGI   = nipa395(5,T)'./nipa394(5,T)'*100;
rSTX  = (nipa32(7,T)+sum(nipa33([9,11],T)))'./nipa119(4,T)'*100;

T     = 69:250; % 1969:1-2014:2
nKCD  = btab100d(T,9)/1000;
nDCD  = atab10d(T,27)/1000;
rKCD  = nKCD./pCD;
rDCD  = nDCD./pCD;

T     = 85:266;  % 1969:1-2014:2
Pop   = (10^3*(civpop(T,2)-civpop(T,3))+10^6*armed(T,2));
H     = hours(T,2)*10^9/4;

Y     = rGDP-rSTX+.04*rKCD+rDCD;
C     = rCND+rCS-(rCND+rCS)./(rCND+rCS+rCD).*rSTX+.04*rKCD+rDCD;
X     = rCD+rGPDI+rGI-rCD./(rCND+rCS+rCD).*rSTX;
% X     = rGPDI;

G     = rGC+rEX-rIM;

%t     = (1969.25:0.25:2014.5)';%quarterly
%t     = (1969:1:2014)';%annual
hpc   = H./Pop;
prd   = Y./H*10^9;
ypc   = Y./Pop*10^9;
xpc   = X./Pop*10^9;
gpc   = G./Pop*10^9;
cgpc  = rGC./Pop*10^9;
cpc0 = C./Pop*10^9;
pcepc = rPCE./Pop*10^9;

t   = (1969.25:0.25:2014.5)';
beighty = find(t==1969.25);

ypc  = ypc(beighty:end);
xpc  = xpc(beighty:end);
hpc  = hpc(beighty:end);
gpc  = gpc(beighty:end);

[~,worktemp.lhpo] = hpfilter(log([ypc hpc xpc gpc]),1600);

% Table II A - observables
yyrelstd     = std(worktemp.lhpo(:,1)); % Case: G=GC+NX
yhrelstd     = std(worktemp.lhpo(:,2))/std(worktemp.lhpo(:,1));
yxrelstd     = std(worktemp.lhpo(:,3))/std(worktemp.lhpo(:,1));
ygtrelstd    = std(worktemp.lhpo(:,4))/std(worktemp.lhpo(:,1));
worktemp.tableIIA1o = [yyrelstd; yhrelstd; yxrelstd; ygtrelstd];
 
yyxcorr       = xcorr(worktemp.lhpo(:,1),worktemp.lhpo(:,1),4,'Coef');
yhxcorr       = xcorr(worktemp.lhpo(:,2),worktemp.lhpo(:,1),4,'Coef');
yxxcorr       = xcorr(worktemp.lhpo(:,3),worktemp.lhpo(:,1),4,'Coef');
ygtxcorr      = xcorr(worktemp.lhpo(:,4),worktemp.lhpo(:,1),4,'Coef');
worktemp.tableIIA2o = [yyxcorr'; yhxcorr'; yxxcorr'; ygtxcorr'];
 
% Table II B - observables
yhxcorr       = xcorr(worktemp.lhpo(:,1),worktemp.lhpo(:,2),4,'Coef'); % Case: G=GC+NX
yxxcorr       = xcorr(worktemp.lhpo(:,1),worktemp.lhpo(:,3),4,'Coef');
ygtxcorr      = xcorr(worktemp.lhpo(:,1),worktemp.lhpo(:,4),4,'Coef');
hxxcorr       = xcorr(worktemp.lhpo(:,2),worktemp.lhpo(:,3),4,'Coef');
hgtxcorr      = xcorr(worktemp.lhpo(:,2),worktemp.lhpo(:,4),4,'Coef');
xgtxcorr      = xcorr(worktemp.lhpo(:,3),worktemp.lhpo(:,4),4,'Coef');
worktemp.tableIIBo = [yhxcorr';yxxcorr';ygtxcorr';hxxcorr';...
    hgtxcorr';xgtxcorr'];
save('worktemp.mat','worktemp','-mat');
