function [x,f,g,code,status]=uncmin(x0,func,param)
%UNCMIN  Unconstrained minimization.
%        [x,f,g,code,status]=uncmin(x0,func,param) solves:
%                          min  f(x,param)                      (*)
%                           x
%        where x0 is the initial guess for the optimal vector x and
%        where  func is a string giving the name of the .m file con-
%        taining code for the function f in (*).
%
%        In addition to returning the optimal x vector, f(x,param), and 
%        g(x)=df(x,param)/dx, UNCMIN returns "code" and "status", where:
%
%        code      =  1   if status(1)=1 or status(j)=1, j=2, 3, and 4
%        code      =  2   if a lower point could not be found
%        code      =  3   if the iteration limit is reached
%        code      =  4   if too many large steps are taken
%
%        status(1) =  1   if norm(g)               < epsa (0 otherwise)
%        status(2) =  1   if norm(g)/(1+f(x))      < gtol (0 otherwise)
%        status(3) =  1   if norm(dx)/(1+norm(x))  < xtol (0 otherwise)
%        status(4) =  1   if df/(1+f(x))           < ftol (0 otherwise)
%        (Note: the user can specify tolerances and other parameters
%         by editing this file.)

%        References:
% 
%        [1] Gill, Philip E., Murray, Walter, and Wright, Margaret H., 
%            PRACTICAL OPTIMIZATION, (New York: Academic Press), 1981.
%        [2] Kahaner, David, Moler, Cleve, and Nash, Stephen, NUMERICAL 
%            METHODS AND SOFTWARE, (Englewood Cliffs, N.J.: Prentice
%            Hall), 
%              

%        Ellen McGrattan,  3-6-89
%        Revised  6-790-

if nargin<3; param=0;end;
%-----------------------------------------------------------------------
%                      USER-DEFINED PARAMETERS
%-----------------------------------------------------------------------

x0=x0(:);                              % we assume that x0 is column 
n=length(x0);                          % vector of length n -- make
                                       % the <func>.m argument  n x 1

%--------------------------default values-------------------------------
epsa=1e-16;                            % tolerance for norm of gradient
epsm=2.22e-16;                         % machine epsilon (=sqrt(2)^2/2-1)
tol=1e-10;
xtol=sqrt(tol);                        % tolerance for change in x
ftol=tol;                              % tolerance for change in f
gtol=tol^(1/3);                        % tolerance for relative gradient
steptl=1e-16;%1e-3*sqrt(epsm);                % tolerance for step length
xsize=max(abs(x0),ones(n,1));          % typical x size
fsize=1.0;                             % typical function size
maxits=10000 ;                         % maximum iterations
iagflg=0;                              % iagflg=1 if analytic gradient 
                                       %   provided (0 otherwise)
iahflg=0;                              % iahflg=1 if analytic hessian
                                       %   provided (0 otherwise)
                                       %   NOTE: analytic hessian not
                                       %   yet coded; set iahflg=0
prt=1;                                 % prt=1 if printing to be done at
                                       %   at intermediate steps

%--------------------------user's values--------------------------------

gtol=1e-7;

%-----------------------------------------------------------------------




%-----------------------------------------------------------------------
%                         INITIALIZATIONS
%-----------------------------------------------------------------------

if min(xsize)<0; xsize=abs(xsize); end;
if fsize<0; fsize=abs(fsize); end;
p=zeros(n,1);
iter=0;
iretcd=-1;
sx=ones(n,1)./xsize;
stpsiz=sqrt(sum((x0.^2).*(sx.^2)));
stepmx=max(1e+3*stpsiz,1e+3);
stepmx=10;
ndigit=fix(-log10(epsm));
rnf=max(10^(-ndigit),epsm);
fddev=diag( max(abs(x0),xsize) );
cddev=rnf^(.33)*fddev;
fddev=sqrt(rnf)*fddev;
%fddev=.00001*fddev;
L=diag(sx);

if iagflg==1;                            % initialize gradient
  eval(['[f0,g0]=',func,'(x0,param);'])
else;
  eval(['f0=',func,'(x0,param);'])
  for i=1:n;
    eval(['g0(i,1)=(',func,'(x0+fddev(:,i),param)-f0)/fddev(i,i);'])
  end;
end;


%-----------------------------------------------------------------------
%                    INITIAL CONVERGENCE CHECK
%-----------------------------------------------------------------------

[code,status,icscmx]=umstop(x0,x0,f0,f0,g0,xtol,ftol,gtol,epsa,xsize, ...
                            fsize,iter,maxits,iretcd,0,0);
if code;
  x=x0;
  f=f0;
  g=g0;
  return;
end;      

%-----------------------------------------------------------------------
%                           MAIN LOOP
%-----------------------------------------------------------------------

while ~code;
  iter=iter+1;

  %
  % solve for search direction, p=- (LL')\g 
  %

  invl=inv(L);
  p=-invl'*invl*g0;

  [x,f,p,mxtake,iretcd]=umlnmin(x0,f0,g0,p,func,param,stepmx,steptl,xsize);

  %  
  % use central differences to compute gradient if a lower point is
  %   not found and forward differences had been used
  %

  if (iretcd==1 & iagflg==0);
    iagflg=-1;
    for i=1:n;
      eval(['g0(i,1)=(',func,'(x0+cddev(:,i),param)-', ...
                        func,'(x0-cddev(:,i),param))/','(2*cddev(i,i));'])
    end;
    invl=inv(L);
    p=-invl'*invl*g0;
    [x,f,p,mxtake,iretcd]=umlnmin(x0,f0,g0,p,func,param,stepmx,steptl,xsize);
  end;

  %
  % compute new gradient vector, g
  %

  if iagflg==-1;
    for i=1:n;
      eval(['g(i,1)=(',func,'(x+cddev(:,i),param)-', ...
                       func,'(x-cddev(:,i),param))/','(2*cddev(i,i));'])
    end;
  elseif iagflg==0;
    for i=1:n;
      eval(['g(i,1)=(',func,'(x+fddev(:,i),param)-f)/fddev(i,i);'])
    end;
  else;
    eval(['[tem,g]=',func,'(x,param);']);
  end;

  %
  % check convergence
  %

  [code,status,icscmx]=umstop(x,x0,f,f0,g,xtol,ftol,gtol,epsa,xsize,fsize, ...
                              iter,maxits,iretcd,mxtake,icscmx);
  %
  % update Hessian = L*L'
  %

  if ~code;

    s=x-x0;
    y=g-g0;
    sy=s'*y;

    if iter==1; 
      L=sqrt(sy/(s'*s))*L; 
      invl=sqrt(s'*s/sy)*invl;
    end;

    if sy>0;
      ls=L'*s;
      tem=L*( eye(n)-(ls*ls')/(ls'*ls)+(invl*y)*(invl*y)'/sy)*L';
      tem=triu(tem)+triu(tem,1)';
      if any(eig(tem) <= -eps) | (norm(tem'-tem,1)/norm(tem,1) > eps) | ...
                                 any(any(imag(tem)))
        if prt; 
          disp('L*S*L not positive definite, eigs are')
          eig(tem)
        end;
      else 
        [Ltem,cholerr] = chol(tem);
        if cholerr~=0;
          if prt; 
            disp('L*S*L not positive definite, eigs are')
            eig(tem)
          end;
        else
          L = Ltem';
        end;
      end;
    end;

    f0=f;
    x0=x;
    g0=g;
  end;

  if prt;    % & rem(iter,50)==0)
    %clc
    %home
    fprintf('Iteration  %g\n',iter)
    disp('   parameter # -- estimate -- gradient')
    disp(sprintf('   %3g %20.8f %20.8f \n',[[1:n]',x,g]'))
    disp(' ')
    disp('   function value (f) -- norm(gradient)/(1+f)')
    disp(sprintf('       %20.8f %20.8f \n',[f,norm(g)/(1+f)]))
    disp(' ')
    disp(' ')
  end
end;
