function [L,Sbar,P0,P,Q,A,B,C,param] = mleq(Theta,extra);
%MLE     Log-likelihood function for the standard growth model with
%        fluctuations in four ``wedges'': z, taul, taux, g.
%        See appendix for details.
%

%        Ellen R. McGrattan, 10-29-02 
%

global ZVAR    % = log([Output, Investment, Hours, Govt Spending ])
               %        (all in per-capita terms)
%---------------------------------------------------------------------
% 1. Default parameters for 
%    gn, gz, beta, delta, psi, sigma, theta, Sbar, P, Q, D, R
%
gn      = (1.015)^(1/4)-1;
gz      = (1.016)^(1/4)-1;
beta    = .9722^(1/4);
delta   = 1-(1-.0464)^(1/4);
psi     = 2.24;
sigma   = 1.000001;
theta   = .35;
adja    = extra(2);
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

ind = [0,0,0,0,0,0,0, 1,1,0,1, ...
       1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1, ...
       1,1,1,1, 0,1,1,1, 0,0,1,1, 0,0,0,1, ...
       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, ...
       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];

Sbar(1) = Theta(1);
Sbar(2) = Theta(2);
Sbar(3) = extra(1);
Sbar(4) = Theta(3);
P(1,1)  = Theta(4);
P(2,1)  = Theta(5);
P(3,1)  = Theta(6);
P(4,1)  = Theta(7);
P(1,2)  = Theta(8);
P(2,2)  = Theta(9);
P(3,2)  = Theta(10);
P(4,2)  = Theta(11);
P(1,3)  = Theta(12);
P(2,3)  = Theta(13);
P(3,3)  = Theta(14);
P(4,3)  = Theta(15);
P(1,4)  = Theta(16);
P(2,4)  = Theta(17);
P(3,4)  = Theta(18);
P(4,4)  = Theta(19);
Q(1,1)  = Theta(20);
Q(2,1)  = Theta(21);
Q(3,1)  = Theta(22);
Q(4,1)  = Theta(23);
Q(2,2)  = Theta(24);
Q(3,2)  = Theta(25);
Q(4,2)  = Theta(26);
Q(3,3)  = Theta(27);
Q(4,3)  = Theta(28);
Q(4,4)  = Theta(29);

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
        1.40, 1.40, 1.40, 1.40, 1.40, 1.40, 1.40, 1.40, 1.40,  ...
              1.40, 1.40, 1.40, 1.40, 1.40, 1.40, 1.40,  ...
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
x     = (eye(4)-P)\P0;
z     = exp(x(1));
taul  = x(2);
taux  = x(3);
g     = exp(x(4));
betah = beta*(1+gz)^(-sigma);
kl    = ((1+taux)*(1-betah*(1-delta))/(betah*theta))^(1/(theta-1))*z;
yk    = (kl/z)^(theta-1);
xi1   = yk-(1+gz)*(1+gn)+1-delta;
xi2   = (1-taul)*(1-theta)*(kl)^theta*z^(1-theta)/psi;
xi3   = xi2/kl;
k     = (xi2+g)/(xi1+xi3);
c     = xi1*k-g;
l     = k/kl;
y     = yk*k;
x     = y-c-g;
r     = theta*y/k;
X0    = [log(k);log(z);taul;taux;log(g);1];
Y0    = [log(y);log(x);log(l);log(g)];

%
% 5b. Approximations for l,y,c,x
%
philh  =-(psi*y*(1-theta)+(1-theta)*(1-taul)*y*(1-l)/l*theta+ ...
         (1-theta)*(1-taul)*y);
philk  = (psi*y*theta+psi*(1-delta)*k-(1-theta)*(1-taul)*y*(1-l)/l *theta)/philh;
philz  = (psi*y*(1-theta)-(1-theta)^2*(1-taul)*y*(1-l)/l)/philh;
phill  = ((1-theta)*(1-taul)*y*(1-l)/l *(1/(1-taul)))/philh;
philx  = 0;
philg  = (-psi*g)/philh;
philkp = (-psi*(1+gz)*(1+gn)*k)/philh;
phiyk  = theta+(1-theta)*philk;
phiyz  = (1-theta)*(1+philz);
phiyl  = (1-theta)*phill;
phiyx  = 0;
phiyg  = (1-theta)*philg;
phiykp = (1-theta)*philkp;
phixk  = -k/x*(1-delta);
phixkp = k/x*(1+gz)*(1+gn);
phixz  = 0;
phixl  = 0;
phixx  = 0;
phixg  = 0;
phick  = y/c*phiyk-x/c*phixk;
phicz  = y/c*phiyz;
phicl  = y/c*phiyl;
phicx  = 0;
phicg  = y/c*phiyg-g/c;
phickp = y/c*phiykp-x/c*phixkp;

%
% 5c. Compute gamma_k
%
coef1  = betah*r*(1-theta);
coef2  = -(1+taux)*psi*(1-sigma)*l/(1-l);
coef3  = -(1+taux)*sigma;
coef4  = betah*(1+taux)*(1+gz)*(1+gn)*adja*adjb;

q      = roots([coef1*philkp+coef2*philkp+coef3*phickp+coef4*phixkp, ...
                coef1*(philk-1)+coef2*(philk-philkp)+coef3*(phick-phickp)+ ...
                coef4*(phixk-1)-(1+taux)*adja*adjb*phixkp, ...
               -coef2*philk-coef3*phick-(1+taux)*adja*adjb*(phixk-1)]);
i      = find(abs(q)<1);
gammak = q(i);

%
% 5d. Compute other gamma's
%
coef  = -coef1-coef2;
kap0  = coef2*[philz;phill;philx;philg]+ ...
        coef3*[phicz;phicl;phicx;phicg]+[0;0;1;0];
kap1  = coef2*philkp+coef3*phickp+coef*(philk+philkp*gammak)- ...
        coef3*(phick+phickp*gammak)+coef1+ ...
        (1+taux)*adja*adjb*phixkp-coef4*(phixk+phixkp*gammak-1);
zet0  = coef*[philz;phill;philx;philg]- ...
        coef3*[phicz;phicl;phicx;phicg]-[coef1;0;betah*(1-delta);0];
zet1  = coef*philkp-coef3*phickp-coef4*phixkp;

gamma = -(kap1*eye(4)+zet1*P')\(kap0+P'*zet0);
gammaz= gamma(1);
gammal= gamma(2);
gammax= gamma(3);
gammag= gamma(4);
gamma0= (1-gammak)*log(k)-gammaz*log(z)-gammal*taul-gammax*taux-gammag*log(g);
gamma = [gammak;gamma;gamma0];
gamhld= gamma;

%
% 5e. State-space system:   X[t+1] = A X[t] + B eps[t+1]
%                           Y[t]   = C X[t] + ome[t]
%                           ome[t] = D ome[t-1] + eta[t],  E eta eta' ~ N(0,R)
%
A     = [ gammak, gammaz, gammal, gammax, gammag, gamma0;
         [0;0;0;0],            P,                     P0;
               0,      0,      0,      0,      0,      1];
B     = [ 0,0,0,0;
            Q;
          0,0,0,0];
C     = [ [phiyk,phiyz,phiyl,phiyx,phiyg]+phiykp*gamma(1:5)';
          [phixk,phixz,phixl,phixx,phixg]+phixkp*gamma(1:5)';
          [philk,philz,phill,philx,philg]+philkp*gamma(1:5)';
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

