%BCA_Simul   Generate time series for a prototype business cycle model 
%            with variations in four wedges. fiscal variables. Before 
%            running this code make sure to first run bca_steady and 
%            bca_wedges. See details in the Appendix to:
%
%               Unmeasured Investment and the 1990s US Hours Boom
%               by Ellen McGrattan and Ed Prescott
                                                                                
%            Ellen McGrattan, 11-1-05
%            Revised, ERM, 3-8-16
%---------------------------------------------------------------------
% 1. Default parameters 
%

bca_params2

T0         = 2008.25;
%bhat       = beta*(1+gamma)^(-sigma);
thet1      = 1-theta;
delt1      = 1-delta;
grate      = (1+eta)*(1+gamma);
adjb       = grate-1+delta;
param      = [adja;beta;delta;eta;gamma;psi;sigma;theta];

%---------------------------------------------------------------------
% 2.  Compute equilibrium
%---------------------------------------------------------------------

load inputs.dat
T          = length(inputs);
% clc
% disp('Business Cycle Accounting for 1-Sector Growth Model')
% disp('---------------------------------------------------')
% iexog = input('  Set any exogenous variables to initial level? (1=yes,0=no) ');
% if iexog==1;
%   disp('    This variable set to initial level? (1=yes,0=no)')
%   labels     = ['    Efficiency wedge?                 ';
%                 '    Labor wedge?                      ';
%                 '    Investment wedge?                 ';
%                 '    Government wedge?                 '];
%   for i=1:4;
%     ifix     = input(labels(i,:));
%     if ifix==1; 
%        inputs(:,i) = inputs(1,i)*ones(T,1); 
%     end;
%   end;
% end;

inputsno = repmat(inputs(1,:),T,1);
inputswz = [inputs(:,1)     inputsno(:,2:4)];
inputswl = [inputsno(:,1)   inputs(:,2) inputsno(:,3:4)];
inputswx = [inputsno(:,1:2) inputs(:,3) inputsno(:,4)];
inputswg = [inputsno(:,1:3) inputs(:,4)];


%%% no wedges sims %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputs = inputsno;
    
o          = ones(T,1);
ewedget    = inputs(:,1);
lwedget    = inputs(:,2);
xwedget    = inputs(:,3);
gwedget    = inputs(:,4);

load vinitial.dat
V0         = vinitial;
deriv      = V0(T,1)/V0(T-1,1);
res        = bca_resid2(V0(:),[param;ks;T;inputs(:);deriv]);
V          = secant('bca_resid2',V0(:),[param;ks;T;inputs(:);deriv]);
V1         = zeros(T,2);
V1(:)      = V;

ktp        = V1(:,1);
ht         = V1(:,2);
kt         = [ks;ktp(1:T-1)];
for i=1:T;
  tem      = roots([adja,-2*(1+adja*adjb),2*(grate*ktp(i)/kt(i)-1+delta)+ ...
             adja*adjb^2]);
  xt(i,1)  = min(tem)*kt(i);
end;
phit       = adja*(xt./kt-adjb).^2/2;
yt         = ewedget.*kt.^theta.*ht.^thet1;
ct         = yt-xt-gwedget;
lt         = 1-ht;
                                                                                
% load observables.dat
% ydat       = observables(:,2)*gdps;
% cdat       = observables(:,4)*gdps;
% hdat       = observables(:,7);
% xdat       = observables(:,5)*gdps;
% kdat(1)    = ks;

load worktemp.mat;
t0 = worktemp.bind;
tT = worktemp.wend;
t             = worktemp.time(t0:tT,1);
gdp           = exp(worktemp.Y(t0:tT,1));%(:,2)*gdps;
ydat          = gdp;
cdat          = exp(worktemp.Y(t0:tT,5));%(:,4)*gdps;
xdat          = exp(worktemp.Y(t0:tT,2));%observables(:,5)*gdps;
hdat          = exp(worktemp.Y(t0:tT,3));%observables(:,7);
kdat(1)       = ks;

for i=1:T-1;
  kdat1(i,1)  = (delt1*kdat(i)+xdat(i))/grate;
  kdat(i+1,1) = kdat1(i);
end;
kdat1(T)   = (delt1*kdat(T)+xdat(T))/grate;
mtfpt      = gdps*ydat./(kdat.^theta.*hdat.^thet1);
prdat      = ydat./hdat;
T1         = T0+T-1;

pftemp.yt= ydat/ydat(1)*100;
pftemp.ht= hdat/hdat(1)*100;
pftemp.xt= xdat/xdat(1)*100;

pftemp.mnoy = yt/yt(1)*100;
pftemp.mnoh = ht/ht(1)*100;
pftemp.mnox = xt/xt(1)*100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% z wedges sims %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputs = inputswz;
    
o          = ones(T,1);
ewedget    = inputs(:,1);
lwedget    = inputs(:,2);
xwedget    = inputs(:,3);
gwedget    = inputs(:,4);

load vinitial.dat
V0         = vinitial;
deriv      = V0(T,1)/V0(T-1,1);
res        = bca_resid2(V0(:),[param;ks;T;inputs(:);deriv]);
V          = secant('bca_resid2',V0(:),[param;ks;T;inputs(:);deriv]);
V1         = zeros(T,2);
V1(:)      = V;

ktp        = V1(:,1);
ht         = V1(:,2);
kt         = [ks;ktp(1:T-1)];
for i=1:T;
  tem      = roots([adja,-2*(1+adja*adjb),2*(grate*ktp(i)/kt(i)-1+delta)+ ...
             adja*adjb^2]);
  xt(i,1)  = min(tem)*kt(i);
end;
phit       = adja*(xt./kt-adjb).^2/2;
yt         = ewedget.*kt.^theta.*ht.^thet1;
ct         = yt-xt-gwedget;
lt         = 1-ht;
                                                                                
% load observables.dat
% ydat       = observables(:,2)*gdps;
% cdat       = observables(:,4)*gdps;
% hdat       = observables(:,7);
% xdat       = observables(:,5)*gdps;
% kdat(1)    = ks;

load worktemp.mat;
t0 = worktemp.bind;
tT = worktemp.wend;
t             = worktemp.time(t0:tT,1);
gdp           = exp(worktemp.Y(t0:tT,1));%(:,2)*gdps;
ydat          = gdp;
cdat          = exp(worktemp.Y(t0:tT,5));%(:,4)*gdps;
xdat          = exp(worktemp.Y(t0:tT,2));%observables(:,5)*gdps;
hdat          = exp(worktemp.Y(t0:tT,3));%observables(:,7);
kdat(1)       = ks;

for i=1:T-1;
  kdat1(i,1)  = (delt1*kdat(i)+xdat(i))/grate;
  kdat(i+1,1) = kdat1(i);
end;
kdat1(T)   = (delt1*kdat(T)+xdat(T))/grate;
mtfpt      = gdps*ydat./(kdat.^theta.*hdat.^thet1);
prdat      = ydat./hdat;
T1         = T0+T-1;

pftemp.mzy = yt/yt(1)*100;
pftemp.mzh = ht/ht(1)*100;
pftemp.mzx = xt/xt(1)*100;
mzye = pftemp.yt - pftemp.mzy;
mzhe = pftemp.ht - pftemp.mzh;
mzxe = pftemp.xt - pftemp.mzx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% taul wedges sims %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputs = inputswl;
    
o          = ones(T,1);
ewedget    = inputs(:,1);
lwedget    = inputs(:,2);
xwedget    = inputs(:,3);
gwedget    = inputs(:,4);

load vinitial.dat
V0         = vinitial;
deriv      = V0(T,1)/V0(T-1,1);
res        = bca_resid2(V0(:),[param;ks;T;inputs(:);deriv]);
V          = secant('bca_resid2',V0(:),[param;ks;T;inputs(:);deriv]);
V1         = zeros(T,2);
V1(:)      = V;

ktp        = V1(:,1);
ht         = V1(:,2);
kt         = [ks;ktp(1:T-1)];
for i=1:T;
  tem      = roots([adja,-2*(1+adja*adjb),2*(grate*ktp(i)/kt(i)-1+delta)+ ...
             adja*adjb^2]);
  xt(i,1)  = min(tem)*kt(i);
end;
phit       = adja*(xt./kt-adjb).^2/2;
yt         = ewedget.*kt.^theta.*ht.^thet1;
ct         = yt-xt-gwedget;
lt         = 1-ht;
                                                                                
% load observables.dat
% ydat       = observables(:,2)*gdps;
% cdat       = observables(:,4)*gdps;
% hdat       = observables(:,7);
% xdat       = observables(:,5)*gdps;
% kdat(1)    = ks;

load worktemp.mat;
t0 = worktemp.bind;
tT = worktemp.wend;
t             = worktemp.time(t0:tT,1);
gdp           = exp(worktemp.Y(t0:tT,1));%(:,2)*gdps;
ydat          = gdp;
cdat          = exp(worktemp.Y(t0:tT,5));%(:,4)*gdps;
xdat          = exp(worktemp.Y(t0:tT,2));%observables(:,5)*gdps;
hdat          = exp(worktemp.Y(t0:tT,3));%observables(:,7);
kdat(1)       = ks;

for i=1:T-1;
  kdat1(i,1)  = (delt1*kdat(i)+xdat(i))/grate;
  kdat(i+1,1) = kdat1(i);
end;
kdat1(T)   = (delt1*kdat(T)+xdat(T))/grate;
mtfpt      = gdps*ydat./(kdat.^theta.*hdat.^thet1);
prdat      = ydat./hdat;
T1         = T0+T-1;

pftemp.mly = yt/yt(1)*100;
pftemp.mlh = ht/ht(1)*100;
pftemp.mlx = xt/xt(1)*100;
mlye = pftemp.yt - pftemp.mly;
mlhe = pftemp.ht - pftemp.mlh;
mlxe = pftemp.xt - pftemp.mlx;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% taux wedges sims %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputs = inputswx;
    
o          = ones(T,1);
ewedget    = inputs(:,1);
lwedget    = inputs(:,2);
xwedget    = inputs(:,3);
gwedget    = inputs(:,4);

load vinitial.dat
V0         = vinitial;
deriv      = V0(T,1)/V0(T-1,1);
res        = bca_resid2(V0(:),[param;ks;T;inputs(:);deriv]);
V          = secant('bca_resid2',V0(:),[param;ks;T;inputs(:);deriv]);
V1         = zeros(T,2);
V1(:)      = V;

ktp        = V1(:,1);
ht         = V1(:,2);
kt         = [ks;ktp(1:T-1)];
for i=1:T;
  tem      = roots([adja,-2*(1+adja*adjb),2*(grate*ktp(i)/kt(i)-1+delta)+ ...
             adja*adjb^2]);
  xt(i,1)  = min(tem)*kt(i);
end;
phit       = adja*(xt./kt-adjb).^2/2;
yt         = ewedget.*kt.^theta.*ht.^thet1;
ct         = yt-xt-gwedget;
lt         = 1-ht;
                                                                                
% load observables.dat
% ydat       = observables(:,2)*gdps;
% cdat       = observables(:,4)*gdps;
% hdat       = observables(:,7);
% xdat       = observables(:,5)*gdps;
% kdat(1)    = ks;

load worktemp.mat;
t0 = worktemp.bind;
tT = worktemp.wend;
t             = worktemp.time(t0:tT,1);
gdp           = exp(worktemp.Y(t0:tT,1));%(:,2)*gdps;
ydat          = gdp;
cdat          = exp(worktemp.Y(t0:tT,5));%(:,4)*gdps;
xdat          = exp(worktemp.Y(t0:tT,2));%observables(:,5)*gdps;
hdat          = exp(worktemp.Y(t0:tT,3));%observables(:,7);
kdat(1)       = ks;

for i=1:T-1;
  kdat1(i,1)  = (delt1*kdat(i)+xdat(i))/grate;
  kdat(i+1,1) = kdat1(i);
end;
kdat1(T)   = (delt1*kdat(T)+xdat(T))/grate;
mtfpt      = gdps*ydat./(kdat.^theta.*hdat.^thet1);
prdat      = ydat./hdat;
T1         = T0+T-1;

pftemp.mxy = yt/yt(1)*100;
pftemp.mxh = ht/ht(1)*100;
pftemp.mxx = xt/xt(1)*100;
mxye = pftemp.yt - pftemp.mxy;
mxhe = pftemp.ht - pftemp.mxh;
mxxe = pftemp.xt - pftemp.mxx;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% gt wedges sims %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputs = inputswg;
    
o          = ones(T,1);
ewedget    = inputs(:,1);
lwedget    = inputs(:,2);
xwedget    = inputs(:,3);
gwedget    = inputs(:,4);

load vinitial.dat
V0         = vinitial;
deriv      = V0(T,1)/V0(T-1,1);
res        = bca_resid2(V0(:),[param;ks;T;inputs(:);deriv]);
V          = secant('bca_resid2',V0(:),[param;ks;T;inputs(:);deriv]);
V1         = zeros(T,2);
V1(:)      = V;

ktp        = V1(:,1);
ht         = V1(:,2);
kt         = [ks;ktp(1:T-1)];
for i=1:T;
  tem      = roots([adja,-2*(1+adja*adjb),2*(grate*ktp(i)/kt(i)-1+delta)+ ...
             adja*adjb^2]);
  xt(i,1)  = min(tem)*kt(i);
end;
phit       = adja*(xt./kt-adjb).^2/2;
yt         = ewedget.*kt.^theta.*ht.^thet1;
ct         = yt-xt-gwedget;
lt         = 1-ht;
                                                                                
% load observables.dat
% ydat       = observables(:,2)*gdps;
% cdat       = observables(:,4)*gdps;
% hdat       = observables(:,7);
% xdat       = observables(:,5)*gdps;
% kdat(1)    = ks;

load worktemp.mat;
t0 = worktemp.bind;
tT = worktemp.wend;
t             = worktemp.time(t0:tT,1);
gdp           = exp(worktemp.Y(t0:tT,1));%(:,2)*gdps;
ydat          = gdp;
cdat          = exp(worktemp.Y(t0:tT,5));%(:,4)*gdps;
xdat          = exp(worktemp.Y(t0:tT,2));%observables(:,5)*gdps;
hdat          = exp(worktemp.Y(t0:tT,3));%observables(:,7);
kdat(1)       = ks;

for i=1:T-1;
  kdat1(i,1)  = (delt1*kdat(i)+xdat(i))/grate;
  kdat(i+1,1) = kdat1(i);
end;
kdat1(T)   = (delt1*kdat(T)+xdat(T))/grate;
mtfpt      = gdps*ydat./(kdat.^theta.*hdat.^thet1);
prdat      = ydat./hdat;
T1         = T0+T-1;

pftemp.mgy = yt/yt(1)*100;
pftemp.mgh = ht/ht(1)*100;
pftemp.mgx = xt/xt(1)*100;
mgye = pftemp.yt - pftemp.mgy;
mghe = pftemp.ht - pftemp.mgh;
mgxe = pftemp.xt - pftemp.mgx;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w1yerrors = [mzye mlye mxye mgye];
temp1 = mean(w1yerrors(:,1:4).^2,1);
pftemp.w1yfz = (1/temp1(1))/sum(temp1.^-1);
pftemp.w1yfl = (1/temp1(2))/sum(temp1.^-1);
pftemp.w1yfx = (1/temp1(3))/sum(temp1.^-1);
pftemp.w1yfg = (1/temp1(4))/sum(temp1.^-1);

w1herrors = [mzhe mlhe mxhe mghe];
temp2 = mean(w1herrors(:,1:4).^2,1);
pftemp.w.w1hfz = (1/temp2(1))/sum(temp2.^-1);
pftemp.w.w1hfl = (1/temp2(2))/sum(temp2.^-1);
pftemp.w.w1hfx = (1/temp2(3))/sum(temp2.^-1);
pftemp.w.w1hfg = (1/temp2(4))/sum(temp2.^-1);

w1xerrors = [mzxe mlxe mxxe mgxe];
temp3 = mean(w1xerrors(:,1:4).^2,1);
pftemp.w.w1xfz = (1/temp3(1))/sum(temp3.^-1);
pftemp.w.w1xfl = (1/temp3(2))/sum(temp3.^-1);
pftemp.w.w1xfx = (1/temp3(3))/sum(temp3.^-1);
pftemp.w.w1xfg = (1/temp3(4))/sum(temp3.^-1);

save pftemp.mat pftemp -mat


% figure(1)
% set(gca,'Fontsize',14)
% plot(T0:T1,hdat/hdat(1)*100,'Linewidth',2,'Color','k')
% line(T0:T1,ht/ht(1)*100,'Linewidth',2,'Color','r','Marker','^')
% ylabel('Index')
% title( 'US Per Capita Hours and Model Prediction')
% legend('US Hours','Model Hours','Location','Southwest')

% figure(2)
% set(gca,'Fontsize',14)
% plot(T0:T1,ydat/ydat(1)*100,'Linewidth',2,'Color','k')
% line(T0:T1,yt/yt(1)*100,'Linewidth',2,'Color','r','Marker','^')
% ylabel('Index')
% % title( 'US Per Capita Real GDP and Model Prediction (/1.02^t)')
% title( 'Spain Per Capita Real GDP and Model Prediction')
% % legend('US GDP','Model GDP','Location','Southwest')
% legend('Spain GDP','Model GDP','Location','Southwest')


% figure(3)
% set(gca,'Fontsize',14)
% plot(T0:T1,xdat/xdat(1)*100,'Linewidth',2,'Color','k')
% line(T0:T1,xt/xt(1)*100,'Linewidth',2,'Color','r','Marker','^')
% ylabel('Index')
% title( 'US Per Capita Real Investment and Model Prediction (/1.02^t)')
% legend('US Investment','Model Investment', ...
%        'Location','Southwest')
% 
% figure(4)
% set(gca,'Fontsize',14)
% plot(T0:T1,cdat/cdat(1)*100,'Linewidth',2,'Color','k')
% line(T0:T1,ct/ct(1)*100,'Linewidth',2,'Color','r','Marker','^')
% ylabel('Index')
% title( 'US Per Capita Real Consumption and Model Prediction (/1.02^t)')
% legend('US Consumption','Model Consumption','Location','Southeast')
% 
% figure(5)
% set(gca,'Fontsize',14)
% plot(T0:T1,hdat/hdat(1)*100,'Linewidth',2,'Color','k')
% line(T0:T1,ewedget/ewedget(1)*100,'Linewidth',2,'Color','r','Marker','^')
% line(T0:T1,lwedget/lwedget(1)*100,'Linewidth',2,'Color','b','Marker','*')
% line(T0:T1,xwedget/xwedget(1)*100,'Linewidth',2,'Color','c','Marker','+')
% title( 'SPA Hours and the BCA Wedges','Fontsize',14)
% legend('SPA Hours','TFP','Labor Wedge','Investment Wedge', ...
%        'Location','Northwest')