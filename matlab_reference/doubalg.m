function [K,S]=double(A,C,Q,R)
%DOUBLE uses the "doubling algorithm" to calculate the Kalman gain, K, and
%  the stationary covariance matrix, S, obtained from the Kalman filter with
%  the following system:
%
%		x[t+1] = Ax[t] + Bu[t] + w1[t+1]
%                 y[t] = Cx[t] + Du[t] + w2[t]
%
%               E [w1(t+1)] [w1(t+1)]' =  [Q  0;
%                 [ w2(t) ] [ w2(t) ]      0  R]
%
%  where x is the mx1 vector of states, u is the nx1 vector of controls, y is
%  the px1 vector of observables, A is mxm, B is mxn, C is pxm, Q is mxm,
%  and R is pxp. 
%
%  INPUTS: A, C, Q, R (not B)
%
%  OUTPUTS: 
%    (i) K= A*S*C'/(R+C*S*C'), the Kalman gain matrix
%   (ii) S= E[x(t)-E[x(t)|y(t-1),y(t-2),..]]*[x(t)-E[x(t)|y(t-1),y(t-2),..]]'
%
m=max(size(A));
v=eye(m);
a0=A';
b=C'*(R\C);
g=Q;
k0=A*g*C'/(C*g*C'+R);
dd=1;
it=1;
maxit=1000;
tol=1e-10;                         % reset tol if necessary
while (dd>tol & it<=maxit);
  v1=(v+b*g);
  vv=v1\a0;
  a1=a0*vv;
  b=b+a0*(v1\(b*a0'));
  g=g+a0'*g*vv;
  k1=A*g*C'/(C*g*C'+R);
  dd=max(max(abs( (k1-k0)./(1+abs(k0)) )));
  it=it+1;
  a0=a1;
  k0=k1;
end
K=k1; S=g;
if it>=maxit; disp('WARNING: Iteration limit of 1000 reached in DOUBLE'); end;
