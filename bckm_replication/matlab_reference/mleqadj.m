function [L,Sbar,P0,P,Q,A,B,C,param] = mleq(Theta,adja);
%MLE     Log-likelihood function for the standard growth model with
%        fluctuations in four ``wedges'': z, taul, taux, g.
%        See appendix for details.
%

%        Ellen R. McGrattan, 10-29-02 
%

load worktemp.mat

param=worktemp.params;
ZVAR  = worktemp.mled(worktemp.iobs:worktemp.eobs,2:5);
% = log([Output, Investment, Hours, Govt Spending ])
               %        (all in per-capita terms)
%---------------------------------------------------------------------
% 1. Default parameters for 
%    gn, gz, beta, delta, psi, sigma, theta, Sbar, P, Q, D, R
%
gn      = param(1);
gz      = param(2);
beta    = param(3);
delta   = param(4);
psi     = param(5);
sigma   = param(6);
theta   = param(7);
adjb    = (1+gn)*(1+gz)-1+delta;
P       = .995*eye(4);
Sbar    = [log(1);.05;.0;log(.07)];
P0      = (eye(4)-P)*Sbar;
Q       = zeros(4);
D       = zeros(4);
R       = zeros(4);
param   = [gn;gz;beta;delta;psi;sigma;theta;adja];

%---------------------------------------------------------------------
% 2. Scale Theta. 

Tsize = ones(length(Theta),1);  
Theta = diag(Tsize)*Theta(:);
%---------------------------------------------------------------------

%---------------------------------------------------------------------
% 3. Assign elements of Theta to underlying parameters.  Use the vector
%    ind to indicate which parameters are to be estimated and which
%    are to be fixed in estimation.  For example, if the vector of
%    parameters is [beta,gamma,delta] and gamma is to be fixed during
%    estimation then set:
%       beta  = Theta(1);
%       delta = Theta(2);
%       ind   = [1,0,1];
%
%    For the current example, the ordering of the parameter vector is:
%       gn, gz, beta, delta, psi, sigma, theta,  Sbar,
%       P(:),  Q(:),  D(:),  R(:)
%        1   2     3      4    5      6      7   8-11
%      12-27  28-43  44-59  60-75

ind = [0,0,0,0,0,0,0, 1,1,1,1, ...
       1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1, ...
       1,1,1,1, 0,1,1,1, 0,0,1,1, 0,0,0,1, ...
       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, ...
       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];

Sbar(1) = Theta(1);
Sbar(2) = Theta(2);
Sbar(3) = Theta(3);
Sbar(4) = Theta(4);
P(1,1)  = Theta(5);
P(2,1)  = Theta(6);
P(3,1)  = Theta(7);
P(4,1)  = Theta(8);
P(1,2)  = Theta(9);
P(2,2)  = Theta(10);
P(3,2)  = Theta(11);
P(4,2)  = Theta(12);
P(1,3)  = Theta(13);
P(2,3)  = Theta(14);
P(3,3)  = Theta(15);
P(4,3)  = Theta(16);
P(1,4)  = Theta(17);
P(2,4)  = Theta(18);
P(3,4)  = Theta(19);
P(4,4)  = Theta(20);
Q(1,1)  = Theta(21);
Q(2,1)  = Theta(22);
Q(3,1)  = Theta(23);
Q(4,1)  = Theta(24);
Q(2,2)  = Theta(25);
Q(3,2)  = Theta(26);
Q(4,2)  = Theta(27);
Q(3,3)  = Theta(28);
Q(4,3)  = Theta(29);
Q(4,4)  = Theta(30);

P0      = (eye(4)-P)*Sbar;

%---------------------------------------------------------------------
% 4. Put bounds on Theta if there are values of the parameters for
%    which L cannot be computed. (Optional)

%    recall ordering: 
%       gn, gz, beta, delta, psi, sigma, theta, Sbar 
%       P(:),  Q(:),  D(:),  R(:)
%        1   2     3      4    5      6      7  8-11
%      12-27  28-43  44-59  60-75
%          
%---------------------------------------------------------------------
Ub  = [.1,.1,.99999,1,3,5,.5,1, 1, 1, 1,  ....
        1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,  ...
              1000, 1000, 1000, 1000, 1000, 1000, 1000,  ...
        1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,  ...
              1000, 1000, 1000, 1000, 1000, 1000, 1000,  ...
        1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,  ...
              1000, 1000, 1000, 1000, 1000, 1000, 1000,  ...
        1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,  ...
              1000, 1000, 1000, 1000, 1000, 1000, 1000]';

Lb  = [0,0,0,0,0,0,0, -1, -1,-1, -5, ...
       -.999,-.999,-.999,-.999,-.999,-.999,-.999,-.999,-.999,  ...
             -.999,-.999,-.999,-.999,-.999,-.999,-.999,  ...
       -1000,-1000,-1000,-1000,-1000,-1000,-1000,-1000,-1000,  ...
             -1000,-1000,-1000,-1000,-1000,-1000,-1000,  ...
       -1000,-1000,-1000,-1000,-1000,-1000,-1000,-1000,-1000,  ...
             -1000,-1000,-1000,-1000,-1000,-1000,-1000,  ...
       -1000,-1000,-1000,-1000,-1000,-1000,-1000,-1000,-1000,  ...
             -1000,-1000,-1000,-1000,-1000,-1000,-1000]';
Ub  = Ub(find(ind));
Lb  = Lb(find(ind));
if (any(Theta>Ub) | any(Theta<Lb));
  L = 1e+20; dL=1e+20*ones(size(Theta));
  return
end;
penalty = 500000*max(max(abs(eig(P)))-.995,0).^2;

%---------------------------------------------------------------------


%---------------------------------------------------------------------
% 5.  Compute equilibrium
%---------------------------------------------------------------------

%
% 5a. Steady state:
%
tem        = (eye(4)-P)\P0;
zs         = exp(tem(1));
tauls      = tem(2);
tauxs      = tem(3);
gs         = exp(tem(4));
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
% 5b. Call subroutine with residuals:
%

Z          = [log(ks);log(ks);log(ks);log(zs);log(zs);tauls;tauls;
              tauxs;tauxs;log(gs);log(gs)];
del        = max(abs(Z)*1e-5,1e-8);
for i=1:11;
  Zp       = Z;
  Zm       = Z;
  Zp(i)    = Z(i)+del(i);
  Zm(i)    = Z(i)-del(i);
  dR(i,1)  = (res_adjust(Zp,param)-res_adjust(Zm,param))/(2*del(i));
end;

%
% 5c. Solution:  log k[t+1] = gamma0 + gammak* log k[t] + gamma* S[t]
%

a0         = dR(1);
a1         = dR(2);
a2         = dR(3);
b0         = dR(4:2:11)';
b1         = dR(5:2:11)';
tem        = roots([a0,a1,a2]);
gammak     = tem(find(abs(tem)<1));
if sum(size(gammak))~=2 || isreal(gammak)==0
    L = 1e+20; dL=1e+20*ones(size(Theta));
    return
else
gamma      = -((a0*gammak+a1)*eye(4)+a0*P')\(b0*P+b1)';
gamma0     = (1-gammak)*log(ks)-gamma'*[log(zs);tauls;tauxs;log(gs)];
Gamma      = [gammak;gamma;gamma0];

%
% 5d. State-space system:   X[t+1] = A X[t] + B eps[t+1]
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
C     = [ [phiyk,phiyz,phiyl,    0,phiyg]+phiykp*Gamma(1:5)';
          [phixk,    0,    0,    0,    0]+phixkp*Gamma(1:5)';
          [philk,philz,phill,    0,philg]+philkp*Gamma(1:5)';
          0,0,0,0,1];
phi0  = Y0-C*X0(1:5);
C     = [C,phi0];
%---------------------------------------------------------------------
% 6. Specify observables (Z). 
%
T          = length(ZVAR);
Y          = log(ZVAR)-log([(1+gz).^[0:T-1]',(1+gz).^[0:T-1]', ...
                             ones(T,1),(1+gz).^[0:T-1]']);
Ybar       = Y(2:T,:)-Y(1:T-1,:)*D';
T          = T-1;
Cbar       = C*A-D*C;
Rbar       = R+C*B*B'*C';
[K,Sigma]  = kfilter(A,Cbar,B*B',Rbar,B*B'*C');
Omega      = Rbar+Cbar*Sigma*Cbar';
Omegai     = inv(Omega);
Xt(1,:)    = X0';
innov(1,:) = Ybar(1,:)-X0'*Cbar';
for i=2:T;
  Xt(i,:)    = Xt(i-1,:)*A'+innov(i-1,:)*K';
  innov(i,:) = Ybar(i,:)-Xt(i,:)*Cbar';
end;
MY         = exp(Xt*Cbar');
DY         = exp(Ybar);
MX         = [exp(Xt(:,1:2)),Xt(:,3:4),exp(Xt(:,5))];
sum1       = innov(1:T,:)'*innov(1:T,:)/T;

L          = 0.5*(T*(log(det(Omega))+trace(Omegai*sum1))+penalty);
end
