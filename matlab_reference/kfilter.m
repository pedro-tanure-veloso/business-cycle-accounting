function [k,s] = kfilter(A,C,V1,V2,V12)

%KFILTER calculates the kalman gain, k, and the stationary covariance 
%  matrix, s, using the kalman filter for:
%  
%		x[t+1] = Ax[t] + Bu[t] + w1[t+1]
%               y[t] = Cx[t] + Du[t] + w2[t]
%
%               E [w1(t+1)] [w1(t+1)]' =  [V1   V12;
%                 [ w2(t) ] [ w2(t) ]      V12' V2 ]
%
%  where x is the mx1 vector of states, u is the nx1 vector of controls, y 
%  is the px1 vector of observables, A is mxm, B is mxn, C is pxm, V1 is 
%  mxm, V2 is pxp, V12 is mxp.
%
%  INPUTS: A, C, V1, V2, V12 (or just A,C,V1,V2 if there are no cross 
%          products, V12=0.
%
%  OUTPUTS:
%    (i) k= (A*s*C'+V12)/(V2+C*s*C'), the Kalman gain matrix
%   (ii) s= A*s*A' + V1 -(A*s*C'+V12)*k' = E[x(t)-E[x(t)|y(t-1),..]]*
%           [x(t)-E[x(t)|y(t-1),..]]'
%
%  NOTES: see also DOUBLE, DRICC, OLRP
%
m=max(size(A));
[rc,cc]=size(C);
if nargin==4; V12=zeros(m,rc); end;
if rank(V2)>=rc;
  A=A-(V12/V2)*C;
  V1=V1-V12*(V2\V12');
  [k,s]=doubalg(A,C,V1,V2);
  k=k+(V12/V2);
else;
  s =zeros(m);
  k0=V12/V2;    %(A*s*C'+V12)/(V2+C*s*C');
  dd=1;
  it=1;
  maxit=1000;
  tol=1e-16;                             % reset tol if necessary Cganhed from -8 to -16
  while (dd>tol & it<=maxit);
    s =A*s*A'+V1-k0*(V2+C*s*C')*k0';
    k1=(A*s*C'+V12)/(V2+C*s*C');
    dd=max(max(abs( (k1-k0)./(1+abs(k0)) )));
    it=it+1;
    k0=k1;
  end;
  k=k1;
  if it>=maxit; 
    disp('WARNING: Iteration limit of 1000 reached in KFILTER.M'); 
  end;
end;

