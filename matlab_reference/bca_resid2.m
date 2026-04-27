function foc=bca_resid(X,param);
%BCA_Resid     Residuals for the prototype business cycle model.
%

%              Ellen R. McGrattan, 11-1-05
%              Revised, ERM, 3-8-16
%______________________________________________________________________________
%
% PARAMETERS

adja       = param(1);
beta       = param(2);
delta      = param(3);
eta        = param(4);
gamma      = param(5);
psi        = param(6);
sigma      = param(7);
theta      = param(8);
k0         = param(9);
T          = param(10);
o          = ones(T,1);
inputs     = zeros(T,4);
inputs(:)  = param(11:10+4*T);
ewedge     = inputs(:,1);
lwedge     = inputs(:,2);
xwedge     = inputs(:,3);
gwedge     = inputs(:,4);
deriv      = param(11+4*T);
bhat       = beta*(1+gamma)^(-sigma);
thet1      = 1-theta;
delt1      = 1-delta;
grate      = (1+eta)*(1+gamma);
adjb       = grate-1+delta;

%______________________________________________________________________________
%
% VARIABLES APPEARING IN RESIDUALS

kn         = X(1:T);
h          = X(T+1:2*T);
k          = [k0;kn(1:T-1)];
for i=1:T;
  tem      = roots([adja,-2*(1+adja*adjb),2*(grate*kn(i)/k(i)-delt1)+ ...
             adja*adjb^2]);
  x(i,1)   = min(tem)*k(i);
end;
phi        = adja*(x./k-adjb).^2/2;
dphi       = adja*(x./k-adjb);
y          = ewedge.*k.^theta.*h.^thet1;
c          = y-x-gwedge;
l          = 1-h;
uc         = c.^(-sigma).*l.^(psi*(1-sigma));
r          = theta*y./k;
t          = 1:T-1;
t1         = 2:T;

foc        = [[uc(t)./(xwedge(t).*(1-dphi(t)))-bhat*uc(t1).*(r(t1)+ ...
               (delt1*o(t1)-phi(t1)+dphi(t1).*x(t1)./k(t1))./ ...
               (xwedge(t1).*(1-dphi(t1))));
                kn(T)-deriv*kn(T-1)];
              psi*c.*h./y-thet1*lwedge.*l]; 
