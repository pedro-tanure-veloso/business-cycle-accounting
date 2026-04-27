%BCA_Wedges  Use parameters from bca_steady and first-order conditions
%            from a prototype business cycle model to back out the  
%            series which are used to form the wedges for business 
%            cycle accounting.  See details in the Appendix to:
%
%               Unmeasured Investment and the 1990s US Hours Boom
%               by Ellen McGrattan and Ed Prescott
                                                                                
%            Ellen McGrattan, 11-1-05
%            Revised, ERM, 3-8-16
                                                                                
%---------------------------------------------------------------------
% Parameters (for steady state parameters, run bca_steady.m first)
%

bca_params2

bhat       = beta*(1+gamma)^(-sigma);
thet1      = 1-theta;
delt1      = 1-delta;
grate      = (1+eta)*(1+gamma);
adjb       = grate-1+delta;

% Load NIPA data 
%
load worktemp.mat
t0 = worktemp.bind;
tT = worktemp.wend;
% t          = observables(:,1);
% gdp        = observables(:,2)*gdps;
% T          = length(gdp);
% y          = gdp;
% c          = observables(:,4)*gdps;
% x          = observables(:,5)*gdps;
% h          = observables(:,7);

t          = worktemp.time(t0:tT,1);
gdp        = exp(worktemp.Y(t0:tT,1));%(:,2)*gdps;
T          = length(gdp);
y          = gdp;
c          = exp(worktemp.Y(t0:tT,5));%(:,4)*gdps;
x          = exp(worktemp.Y(t0:tT,2));%observables(:,5)*gdps;
h          = exp(worktemp.Y(t0:tT,3));%observables(:,7);

%
% Back out variables
%
gwedge     = exp(worktemp.Y(t0:tT,4));%gdp-x-c;
k          = ks*ones(T,1);
kn         = k;
for i=1:T-1;
  phi(i)   = adja*(x(i)/k(i)-adjb)^2/2;
  kn(i)    = (delt1*k(i)+x(i)-phi(i)*k(i))/grate;
  k(i+1)   = kn(i);
end;
phi(T)     = adja*(x(T)/k(T)-adjb)^2/2;
kn(T)      = (delt1*k(T)+x(T)-phi(T)*k(T))/grate;
ewedge     = y./(k.^theta.*h.^thet1);
lwedge     = psi*c.*h./(thet1*y.*(1-h));    
uc         = c.^(-sigma).*(1-h).^(psi*(1-sigma)); 
r          = theta*y./k;
dphi       = adja*(x./k-adjb);
xwedge     = xwedges*ones(T,1);
for i=1:T-1;
  dtem        = (delt1-phi(i+1)+dphi(i+1).*x(i+1)/k(i+1))/(1-dphi(i+1));
  xwedge(i+1) = (bhat*uc(i+1)*dtem)/ ...
                (uc(i)/(xwedge(i).*(1-dphi(i)))-bhat*uc(i+1)*r(i+1));
end;

inputs     = [ewedge,lwedge,xwedge,gwedge];
initial    = [kn,h];

save inputs.dat inputs  -ascii -double
save vinitial.dat initial  -ascii -double

