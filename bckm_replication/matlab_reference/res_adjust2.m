function R = res_wedge2(Z,param);
%RES_WEDGE     residual for the growth model in ``Business Cycle Accounting.''
%

%              Ellen R. McGrattan, 11-1-02
%              Revised, ERM, 11-1-02

%______________________________________________________________________________
%
% PARAMETERS

gn         = param(1);
gz         = param(2);
beta       = param(3);
delta      = param(4);
psi        = param(5);
sigma      = param(6);
theta      = param(7);
a          = param(8);
logzbar    = param(9);
taulbar    = param(10);
tauxbar    = param(11);
loggbar    = param(12);
Az         = param(13);
Al         = param(14);
Ax         = param(15);
Ag         = param(16);
b          = (1+gz)*(1+gn)-1+delta;
beth       = beta*(1+gz)^(-sigma);
                                                                                
%______________________________________________________________________________
%
% VARIABLES APPEARING IN RESIDUALS
                                                                                
k2         = exp(Z(1));
k1         = exp(Z(2));
k          = exp(Z(3));
sz1        = exp(Z(4));
sz         = exp(Z(5));
staul1     = Z(6);
staul      = Z(7);
staux1     = Z(8);
staux      = Z(9);
sg1        = exp(Z(10));
sg         = exp(Z(11));
z1         = exp(Az*log(sz1)+(1-Az)*logzbar);
z          = exp(Az*log(sz)+(1-Az)*logzbar);
taul1      = Al*staul1+(1-Al)*taulbar;
taul       = Al*staul+(1-Al)*taulbar;
taux1      = Ax*staux1+(1-Ax)*tauxbar;
taux       = Ax*staux+(1-Ax)*tauxbar;
g1         = exp(Ag*log(sg1)+(1-Ag)*loggbar);
g          = exp(Ag*log(sg)+(1-Ag)*loggbar);

tem        = roots([a,-2*(1+a*b),2*((1+gz)*(1+gn)*k1/k-1+delta)+a*b^2]);
x          = min(tem)*k;
tem        = roots([a,-2*(1+a*b),2*((1+gz)*(1+gn)*k2/k1-1+delta)+a*b^2]);
x1         = min(tem)*k1;
phi        = a/2*(x/k-b)^2;
phi1       = a/2*(x1/k1-b)^2;
dphi       = a*(x/k-b);
dphi1      = a*(x1/k1-b);
l          = 1/(1+.75*psi/(1-taul)/(1-theta));
l1         = 1/(1+.75*psi/(1-taul1)/(1-theta));
l=.99;
l1=.99;
for i=1:5;
  y        = k^theta*(z*l)^(1-theta);
  y1       = k1^theta*(z1*l1)^(1-theta);
  c        = y-x-g;
  c1       = y1-x1-g1;
  res      = psi*c*l/y-(1-taul)*(1-theta)*(1-l); 
  res1     = psi*c1*l1/y1-(1-taul1)*(1-theta)*(1-l1); 

  lp       = l+.0001;
  l1p      = l1+.0001;
  y        = k^theta*(z*lp)^(1-theta);
  y1       = k1^theta*(z1*l1p)^(1-theta);
  c        = y-x-g;
  c1       = y1-x1-g1;
  dres     = (psi*c*lp/y-(1-taul)*(1-theta)*(1-lp)-res)/.0001; 
  dres1    = (psi*c1*l1p/y1-(1-taul1)*(1-theta)*(1-l1p)-res1)/.0001; 
 
  l        = l-res/dres; 
  l1       = l1-res1/dres1; 
end;
y          = k^theta*(z*l)^(1-theta);
y1         = k1^theta*(z1*l1)^(1-theta);
c          = y-x-g;
c1         = y1-x1-g1;

R          = (1+taux)*c^(-sigma)*(1-l)^(psi*(1-sigma))/(1-dphi)- ...
             beth*c1^(-sigma)*(1-l1)^(psi*(1-sigma))* ...
             (theta*y1/k1+(1-delta-phi1+dphi1*x1/k1)*(1+taux1)/(1-dphi1));

