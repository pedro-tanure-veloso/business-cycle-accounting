function x = secant(func,x0,param,crit,maxit);
%SECANT   Finds the fixed point of a system of equations.
%         x=secant(func,x0,param,crit,maxit)  uses the secant method
%         to find the zeros of the following system:
%
%                        f1(z1,z2,...zn)=0
%                        f2(z1,z2,...zn)=0                (*)
%                          :           :
%                        fn(z1,z2,...zn)=0
%
%         where x=[z1,z2,...,zn]' are those values of zi that solve (*).
%         To use secant, you must specify "func" which is a string giving
%         the name of an m-file.  Do not use "f.m","x.m","x0.m".
%         This m-file takes inputs x and param and returns f, where
%         f = [f1(x,param),f2(x,param),...fn(x,param)]'.  The arguments 
%         x0, crit, and maxit are the initial guess for the fixed point,
%         the convergence criterion, and the iteration limit.

%         Ellen R. McGrattan, 11-13-87
%
if nargin < 4; crit=1e-8; maxit=100; end;
n   = length(x0);
del = diag(max(abs(x0)*1e-4,1e-8));
for i=1:maxit;
  eval(['f        =',func,'(x0,param);']);
  for j=1:n;
    eval(['fj(:,j)  = (f-',func,'(x0-del(:,j),param))/del(j,j);']);
  end;
  x=x0-fj\f;
  if max( abs(x-x0)./(1+abs(x0)) )<crit; 
    break; 
  end;
  x0=x;
end;
if i>=maxit; 
  sprintf('WARNING: iteration limit of %g exceeded for NEWTON',maxit)
end;
