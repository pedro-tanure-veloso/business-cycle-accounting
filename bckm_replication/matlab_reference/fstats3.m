%FSTATS    This code computes the following statistic:
%
%                        1/sum((dly-dly(i))^2)
%              f(i) = --------------------
%                     sum_j (1/(sum(dly-dly(j)^2))
%
%          where dly is logged and filtered data and dly(i) is
%          logged and filtered component of that data due to wedge
%          i=(A,taul,taux,g).
%

%          Ellen McGrattan, 1-9-16

load worktemp.mat;%ADDED BY PEDRO
t = worktemp.time;%ADDED BY PEDRO

ilog        = 0;  % 1--> use logs,  0--> use levels
ifilt       = 0;  % 1--> HP filter, 0--> no filter
i1          = find(t==2008.25); %ADDED BY PEDRO
i2          = find(t==2011.75); %ADDED BY PEDRO
% i1          = find(t==2008.125);
% i2          = find(t==2011.875);
%i1          = find(t==1948.125);
%i2          = find(t==2015.625);

%
% Data and components due to different wedges:
%
% Output      = exp(lyt(i1:i2)-lyt(Y0));
% OC          = [exp(YMz(i1:i2,1)-YMz(Y0,1)),exp(YMl(i1:i2,1)-YMl(Y0,1)), ...
%                exp(YMx(i1:i2,1)-YMx(Y0,1)),exp(YMg(i1:i2,1)-YMg(Y0,1))];
Output      = worktemp.w.yt(i1:i2);
OC          = [worktemp.w.mzy(i1:i2),worktemp.w.mly(i1:i2),...
               worktemp.w.mxy(i1:i2),worktemp.w.mgy(i1:i2)]/100;  %ADDED BY PEDRO

% Labor       = exp(llt(i1:i2)-llt(Y0));
% LC          = [exp(YMz(i1:i2,3)-YMz(Y0,3)),exp(YMl(i1:i2,3)-YMl(Y0,3)), ...
%                exp(YMx(i1:i2,3)-YMx(Y0,3)),exp(YMg(i1:i2,3)-YMg(Y0,3))];
Labor       = worktemp.w.ht(i1:i2);
LC          = [worktemp.w.mzh(i1:i2),worktemp.w.mlh(i1:i2),...
               worktemp.w.mxh(i1:i2),worktemp.w.mgh(i1:i2)]/100;  %ADDED BY PEDRO

% Investment  = exp(lxt(i1:i2)-lxt(Y0));
% XC          = [exp(YMz(i1:i2,2)-YMz(Y0,2)),exp(YMl(i1:i2,2)-YMl(Y0,2)), ...
%                exp(YMx(i1:i2,2)-YMx(Y0,2)),exp(YMg(i1:i2,2)-YMg(Y0,2))];

Investment  = worktemp.w.xt(i1:i2);
XC          = [worktemp.w.mzx(i1:i2),worktemp.w.mlx(i1:i2),...
               worktemp.w.mxx(i1:i2),worktemp.w.mgx(i1:i2)]/100;  %ADDED BY PEDRO

%
% F-statistics
%
if ilog==0;
  dly       = Output;
  dll       = Labor;
  dlx       = Investment;
  dlyc      = OC;
  dllc      = LC;
  dlxc      = XC;
else 
  dly       = log(Output);
  dll       = log(Labor);
  dlx       = log(Investment);
  dlyc      = log(OC);
  dllc      = log(LC);
  dlxc      = log(XC);
end;
if ifilt==1;
  dly       = dly-hptrend(dly,1600);
  dll       = dll-hptrend(dll,1600);
  dlx       = dlx-hptrend(dlx,1600);
  for i=1:4;
    dlyc(:,i)  = dlyc(:,i)-hptrend(dlyc(:,i),1600);
    dllc(:,i)  = dllc(:,i)-hptrend(dllc(:,i),1600);
    dlxc(:,i)  = dlxc(:,i)-hptrend(dlxc(:,i),1600);
  end;
end;

for i=1:4;
 f(1,i)     = 1/sum((dly-dlyc(:,i)+1e-10).^2);
 f(2,i)     = 1/sum((dll-dllc(:,i)+1e-10).^2);
 f(3,i)     = 1/sum((dlx-dlxc(:,i)+1e-10).^2);
end;
fstat       = f./(sum(f')'*ones(1,4));

%
% Print out results
%

disp('                     F-Statistics ')
disp('-------------------------------------------------------')
disp('                    Contribution of wedge:             ')
disp('              Efficiency  Labor   Investment  G.Cons. ')
disp('-------------------------------------------------------')
disp(sprintf(' Output     %9.2f %9.2f %9.2f %9.2f',fstat(1,:)))
disp(sprintf(' Labor      %9.2f %9.2f %9.2f %9.2f',fstat(2,:)))
disp(sprintf(' Investment %9.2f %9.2f %9.2f %9.2f',fstat(3,:)))
disp('-------------------------------------------------------')

