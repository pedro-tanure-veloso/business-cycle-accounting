function x0 = initmle(Sbar)

load worktemp.mat

param=worktemp.params;
ZVAR  = worktemp.Y(worktemp.iobs:worktemp.eobs,1:4);
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
P       = eye(4)*0.995;
P0      = (eye(4)-P)*Sbar;
% Sbar    = [log(1);.05;.0;log(.07)];
% P0      = (eye(4)-P)*Sbar;
% Q       = zeros(4);
% D       = zeros(4);
% R       = zeros(4);

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

% compute sample means
Ym = mean(exp(ZVAR));

x0 = [ys-Ym(1);xs/ys-Ym(2);ls-Ym(3);gs/ys-Ym(4)];

