function [A,B,C,X0] = fixexpadj(Sbar,P,Q,s0,As,params,adja);
%FIXEXP  Compute decision functions holding fixed expectations and
%        setting some of the wedges to constants.
%        See appendix for details.
%

%        Ellen R. McGrattan, 10-6-06 
%

%---------------------------------------------------------------------
% 1. Default parameters for 
%    gn, gz, beta, delta, psi, sigma, theta, Sbar, P, Q, D, R
%
gn      = params(1);
gz      = params(2);
beta    = params(3);
delta   = params(4);
psi     = params(5);
sigma   = params(6);
theta   = params(7);
P0      = (eye(4)-P)*Sbar;

param   = [gn;gz;beta;delta;psi;sigma;theta;adja;s0;As];

%---------------------------------------------------------------------
% 2.  Compute equilibrium
%---------------------------------------------------------------------

%
% 2a. Steady state:
%
zs         = exp(Sbar(1));
tauls      = Sbar(2);
tauxs      = Sbar(3);
gs         = exp(Sbar(4));
beth       = beta*(1+gz)^(-sigma);
kls        = ((1+tauxs)*(1-beth*(1-delta))/(beth*theta))^(1/(theta-1))*zs;
A          = (zs/kls)^(1-theta)-(1+gz)*(1+gn)+1-delta;
B          = (1-tauls)*(1-theta)*kls^theta*zs^(1-theta)/psi;
ks         = (B+gs)/(A+B/kls);
cs         = A*ks-gs;
ls         = ks/kls;
ys         = ks^theta*(zs*ls)^(1-theta);
xs         = ys-cs-gs;
X0         = [log(ks);log(zs);tauls;tauxs;log(gs);1];
Y0         = [log(ys);log(xs);log(ls);log(gs)];

%
% 2b. Call subroutine with residuals:
%

Z          = [log(ks);log(ks);log(ks);log(zs);log(zs);tauls;tauls;
              tauxs;tauxs;log(gs);log(gs)];
del        = max(abs(Z)*1e-5,1e-8);
for i=1:11;
  Zp       = Z;
  Zm       = Z;
  Zp(i)    = Z(i)+del(i);
  Zm(i)    = Z(i)-del(i);
  dR(i,1)  = (res_adjust2(Zp,param)-res_adjust2(Zm,param))/(2*del(i));
end;

%
% 2c. Solution:  log k[t+1] = gamma0 + gammak* log k[t] + gamma* S[t]
%

a0         = dR(1);
a1         = dR(2);
a2         = dR(3);
b0         = dR(4:2:11)';
b1         = dR(5:2:11)';
tem        = roots([a0,a1,a2]);
gammak     = tem(find(abs(tem)<1));
gamma      = -((a0*gammak+a1)*eye(4)+a0*P')\(b0*P+b1)';
gamma0     = (1-gammak)*log(ks)-gamma'*[log(zs);tauls;tauxs;log(gs)];
Gamma      = [gammak;gamma;gamma0];

%
% 2d. State-space system:   X[t+1] = A X[t] + B eps[t+1]
%                           Y[t]   = C X[t] + ome[t]
%                           ome[t] = D ome[t-1] + eta[t],  E eta eta' ~ N(0,R)
%
philh      =-(psi*ys*(1-theta)+(1-theta)*(1-tauls)*ys*(1-ls)/ls*theta+ ...
             (1-theta)*(1-tauls)*ys);
philk      = (psi*ys*theta+psi*(1-delta)*ks- ...
             (1-theta)*(1-tauls)*ys*(1-ls)/ls *theta)/philh;
philz      = (psi*ys*(1-theta)-(1-theta)^2*(1-tauls)*ys*(1-ls)/ls)/philh;
phill      = ((1-theta)*(1-tauls)*ys*(1-ls)/ls *(1/(1-tauls)))/philh;
philg      = (-psi*gs)/philh;
philkp     = (-psi*(1+gz)*(1+gn)*ks)/philh;
phiyk      = theta+(1-theta)*philk;
phiyz      = (1-theta)*(1+philz);
phiyl      = (1-theta)*phill;
phiyg      = (1-theta)*philg;
phiykp     = (1-theta)*philkp;
phixk      = -ks/xs*(1-delta);
phixkp     = ks/xs*(1+gz)*(1+gn);

A     = [ gammak,                         gamma', gamma0;
         [0;0;0;0],            P,                     P0;
               0,      0,      0,      0,      0,      1];
B     = [ 0,0,0,0;
            Q;
          0,0,0,0];
C     = [ [phiyk,phiyz*As(1),phiyl*As(2),    0,phiyg*As(4)]+phiykp*Gamma(1:5)';
          [phixk,       0,       0,    0,                0]+phixkp*Gamma(1:5)';
          [philk,philz*As(1),phill*As(2),    0,philg*As(4)]+philkp*Gamma(1:5)';
          0,0,0,0,As(4)];
phi0  = Y0-C*X0(1:5);
C     = [C,phi0];

