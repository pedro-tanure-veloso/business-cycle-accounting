function [sbarpq,Lk] = runmleadj()

%RUNMLE     Compute MLE estimates for the benchmark model using the
%           sample period 1959:1--2004:3 and bootstrap standard errors 
%           if constraints bind on matrix P. The likelihood function
%           is in file mleq.m
%
                                                                                
%           Ellen McGrattan, 2-16-02
%           Revised, ERM, 3-18-05

load worktemp.mat

sbari = fsolve(@initmle,[0;0.05;0.0;log(0.2)]);

x0a=[ ...
   sbari(1)
   sbari(2)
   sbari(3)
   sbari(4)
   0.995
   0
   0
   0
   0
   0.995
   0
   0
   0
   0
   0.995
   0
   0
   0
   0.
   0.995
   0.01164416023534
   0.00123692627331
  -0.00533935134032
  -0.00079381859001
   0.00644997326129
   0.00230641582338
   0.00639266420784
   0.00981139512775
   0.01374609848920
   0.00560226065257]; % result from runmle (with no adjustment costs)

x0b =[ ...
   sbari(1)
   sbari(2)
   sbari(3)
   sbari(4)
   0.995
   0
   0
   0
   0
   0.995
   0
   0
   0
   0
   0.995
   0
   0
   0
   0.
   0.995
   0.02382866874935
  -0.01028081000637
   0.00553372214162
                  0
   0.02739908940082
  -0.03354291115539
                  0
   0.04951600718821
                  0
   0.10194597840286]; % result from initpw with annual adja=3.22
             
x0c=[ ...
   sbari(1)
   sbari(2)
   sbari(3)
   sbari(4)
   0.995
   0
   0
   0
   0
   0.995
   0
   0
   0
   0
   0.995
   0
   0
   0
   0.
   0.995
   0.02396761427982
  -0.00987436176711
  -0.01693235174207
                  0
   0.02737005516313
  -0.06560608935313
                  0
   0.12084347484485
                  0
   0.10034489721325];  % result from initpw with annual adja=12.88
                       % and maximum eigenvalue on P^(1/4) scaled 
XXX = [x0a x0b x0c];
adjcs = [0 12.88 4*12.88];
                       
x0      = XXX(:,worktemp.adjc);
adja    = adjcs(worktemp.adjc);

pb = worktemp.optimnum.pb;
nps = worktemp.optimnum.nps; %number of periods to check improvement (added)
X       = zeros(30,nps);hh = waitbar(0,'Maximizing Likelihood...0%');
[x1,f,g,code,status]   = uncmin(x0,'mleqadj',adja);
X(:,1) = x1;close(hh);
hh = waitbar(1/nps,['Maximizing Likelihood...',num2str(1/nps*100),'%']);
F(1,1) = f;
[x1,f,g,code,status]   = uncmin(x1,'mleqadj',adja);
X(:,2) = x1;close(hh);
hh = waitbar(2/nps,['Maximizing Likelihood...',num2str(2/nps*100),'%']);
F(2,1) = f;
%
% Move away to see if we get more improvement
%
x2 = x1;

for i=3:nps; %originally 50
  [x2,f,g,code,status]   = uncmin(x2*pb,'mleqadj',adja);
  F(i,1) = f;
  X(:,i) = x2;close(hh);
  hh = waitbar(i/nps,['Maximizing Likelihood...',num2str(i/nps*100),'%']);
  disp(F)
  pause(2)
end;
close(hh)
i           = find(F==min(F));
if sum(size(i))>2
    i = i(1);
end
x           = X(:,i);
L           = F(i);
% ibind       = 1;
param       = zeros(30,1);
ind         = 1:30;
param(ind)  = x;

worktemp.mlemax.F = F;
worktemp.mlemax.X = X;
worktemp.mlemax.X0 = x0;
worktemp.mlemax.nps = nps;
worktemp.mlemax.optind = i;
save('worktemp.mat','worktemp','-mat');

% if ibind==0;
%   %
%   % Standard errors for the case with nonbinding constraints
%   %
%   del   = diag(max(abs(param)*1e-4,1e-8));
%   for i=1:length(param);
%     [f1,f2]  = mleseqadj(param+del(:,i),adja);
%     [m1,m2]  = mleseqadj(param-del(:,i),adja);
%     dL(i,1)  = (f1-m1)/(2*del(i,i));
%     dLt(i,:) = (f2-m2)'/(2*del(i,i));
%   end;
%   [n,m] = size(dLt);
%   sum1  = 0;
%   for t=1:m;
%     sum1=sum1+dLt(:,t)*dLt(:,t)';
%   end;
%   se = diag(sqrt(inv(sum1)));
% else;
%   %
%   % Bootstrapped standard errors for the case with binding constraints
%   %
%   [L,Lt,ut,X0,Cbar,A,K,D,gz] = mleseqadj(param,adja);
%   B         = 500;   % number of bootstrap replications
%   T         = length(ut);
%   Ybar      = 0*ut;
%   nx        = length(X0);
%   Xt        = zeros(T+1,nx);
%   Y         = zeros(T+1,4);
%   Y(1,:)    = log(ZVAR(1,:));
%   Data      = zeros(4*(T+1),B);
%   Theta     = zeros(B,length(ind));
%   thet0     = param(ind);
%   rand('state',100462)
%   for i=1:B;
%     %
%     % Draw u's uniformly from sample {ut(1),ut(2)....ut(T)}
%     % 
%     Xt(1,:)   = X0';
%     for j=1:T;
%       k             = ceil(rand*T);  % draw number between 1 and T
%       Ybar(j,:)     = Xt(j,:)*Cbar'+ut(k,:);
%       Xt(j+1,:)     = Xt(j,:)*A'+ut(k,:)*K';
%     end;
%     for j=2:T+1;
%       Y(j,:) = Y(j-1,:)*D'+Ybar(j-1,:);
%     end;
%     ZVAR      = exp(Y+log([(1+gz).^[0:T]',(1+gz).^[0:T]', ...
%                              ones(T+1,1),(1+gz).^[0:T]']));
%     Data(:,i) = ZVAR(:);
%     theta     = uncmin(thet0,'mleqadj',adja);
%     if theta(21)<0;
%       theta(21:24) = -theta(21:24);
%     end;
%     if theta(25)<0;
%       theta(25:27) = -theta(25:27);
%     end;
%     if theta(28)<0;
%       theta(28:29) = -theta(28:29);
%     end;
%     if theta(30)<0;
%       theta(30)    = -theta(30);
%     end;
%     Theta(i,:) = theta';
%   end;
%   se  = std(Theta)';
% end;
%
% Print results 
%

disp('Results from Maximum Likelihood Estimation')
disp('------------------------------------------')
disp(' ')
% disp('  [Theta, Standard Errors] ')
disp('  [Theta] ')
% disp(sprintf(' %10.3e %10.3e\n', [param(ind),se]'))
disp(sprintf(' %10.3e \n', [param(ind)]'))
disp(' ')
fprintf('  L(Theta) = %g ',L)
disp(' ')
                                                                                
sbarpq=param(ind); Lk=L;

