function [] = gwedges2()
% This function takes as input a 3 letter country code and backs out the
% wedges given the results of the MLE estimation.

% Assumes data folder to be current folder.

% resets path
%restoredefaultpath;

% loads country input structure file
load('worktemp.mat');
mle = worktemp.mle;
qadj  = 0.25/sum(worktemp.params([1:2 4]));
adjcs = [0 qadj 4*qadj];
adja    = adjcs(worktemp.adjc);
t    = worktemp.time;
ZVAR = mle.obs;
x0   = mle.Theta;
   
[L,Sbar,P0,P,Q,A,B,C,param] = mleqadj(x0,adja);
Y0       = worktemp.bind;
gn       = param(1);
gz       = param(2);
beta     = param(3);
delta    = param(4);
psi      = param(5);
sigma    = param(6);
theta    = param(7);
z        = exp(Sbar(1));
taul     = Sbar(2);
taux     = Sbar(3);
g        = exp(Sbar(4));
betah    = beta*(1+gz)^(-sigma);
kl       = ((1+taux)*(1-betah*(1-delta))/(betah*theta))^(1/(theta-1))*z;
yk       = (kl/z)^(theta-1);
xi1      = yk-(1+gz)*(1+gn)+1-delta;
xi2      = (1-taul)*(1-theta)*(kl)^theta*z^(1-theta)/psi;
k        = (xi2+g)/(xi1+xi2/kl);
c        = xi1*k-g;
l        = k/kl;
y        = yk*k;
x        = y-c-g;
lk       = log(k);
lc       = log(c);
ll       = log(l);
ly       = log(y);
lx       = log(x);
lg       = log(g);
lz       = log(z);
T        = length(ZVAR);
Y        = log(ZVAR)-log([(1+gz).^[0:T-1]',(1+gz).^[0:T-1]', ...
                         ones(T,1),(1+gz).^[0:T-1]',(1+gz).^[0:T-1]']);
lyt      = Y(:,1);
lxt      = Y(:,2);
llt      = Y(:,3);
lgt      = Y(:,4);
lgct     = Y(:,5);
% lnxt     = Y(:,6);
lkt(1,1) = lk;
Kt(1,1)  = exp(lk);

for i=1:T;
  lktp(i,1)  = lk+((1-delta)*(lkt(i)-lk)+x/k*(lxt(i)-lx))/(1+gz)/(1+gn);
  lkt(i+1,1) = lktp(i);
  Ktp(i,1)   = ((1-delta)*Kt(i)+exp(lxt(i)))/(1+gz)/(1+gn);
  Kt(i+1,1)  = Ktp(i);
end;
lkt      = lkt(1:T);
Kt       = Kt(1:T);
lct      = lc+(y*(lyt-ly)-x*(lxt-lx)-g*(lgt-lg))/c;
lzt      = lz+(lyt-ly-theta*(lkt-lk))/(1-theta)-llt+ll;
tault    = taul+(1-taul)*(lyt-ly-lct+lc-1/(1-l)*(llt-ll));
tauxt    = (lxt-C(2,1)*lkt-C(2,2)*lzt-C(2,3)*tault-C(2,5)*lgt-C(2,6))/C(2,4);
tauxchk  = (lyt-C(1,1)*lkt-C(1,2)*lzt-C(1,3)*tault-C(1,5)*lgt-C(1,6))/C(1,4);
Ct       = exp(lyt)-exp(lxt)-exp(lgt);
Zt       = (exp(lyt)./(Kt.^theta.*exp(llt).^(1-theta))).^(1/(1-theta));
Tault    = 1-psi/(1-theta)* (Ct./exp(lyt)) .*(exp(llt)./(1-exp(llt)));

% simulation part
Xt0      = [lkt,lzt,tault,tauxt,lgt,ones(T,1)];
Xt0_gc   = [lkt,lzt,tault,tauxt,lgct,ones(T,1)];
% Xt0_nx   = [lkt,lzt,tault,tauxt,lnxt,ones(T,1)];
YM0      = Xt0*C';
YM0_gc   = Xt0_gc*C';

s0       = [lzt(Y0);tault(Y0);tauxt(Y0);lgt(Y0)];
s0_gc    = [lzt(Y0);tault(Y0);tauxt(Y0);lgct(Y0)];
% s0_nx    = [lzt(Y0);tault(Y0);tauxt(Y0);lnxt(Y0)];

[i,j,C0] = fixexpadj(Sbar,P,Q,s0,[0;0;0;0],param,adja);
[i,j,C1] = fixexpadj(Sbar,P,Q,s0,[1;0;0;0],param,adja);
[i,j,C2] = fixexpadj(Sbar,P,Q,s0,[0;1;0;0],param,adja);
[i,j,C3] = fixexpadj(Sbar,P,Q,s0,[0;0;1;0],param,adja);
[i,j,C4] = fixexpadj(Sbar,P,Q,s0,[0;0;0;1],param,adja);

[igc,jgc,C0gc] = fixexpadj(Sbar,P,Q,s0_gc,[0;0;0;0],param,adja);
[igc,jgc,C1gc] = fixexpadj(Sbar,P,Q,s0_gc,[1;0;0;0],param,adja);
[igc,jgc,C2gc] = fixexpadj(Sbar,P,Q,s0_gc,[0;1;0;0],param,adja);
[igc,jgc,C3gc] = fixexpadj(Sbar,P,Q,s0_gc,[0;0;1;0],param,adja);
[igc,jgc,C4gc] = fixexpadj(Sbar,P,Q,s0_gc,[0;0;0;1],param,adja);

% 
% [inx,jnx,C0nx] = fixexpadj(Sbar,P,Q,s0_nx,[0;0;0;0],param,adja);
% [inx,jnx,C1nx] = fixexpadj(Sbar,P,Q,s0_nx,[1;0;0;0],param,adja);
% [inx,jnx,C2nx] = fixexpadj(Sbar,P,Q,s0_nx,[0;1;0;0],param,adja);
% [inx,jnx,C3nx] = fixexpadj(Sbar,P,Q,s0_nx,[0;0;1;0],param,adja);
% [inx,jnx,C4nx] = fixexpadj(Sbar,P,Q,s0_nx,[0;0;0;1],param,adja);

o        = ones(T,1);

YMn      = (Xt0-o*Xt0(Y0,:))*C0'+o*YM0(Y0,:);
YMz      = (Xt0-o*Xt0(Y0,:))*(C1-C0)'+o*YM0(Y0,:);
YMl      = (Xt0-o*Xt0(Y0,:))*(C2-C0)'+o*YM0(Y0,:);
YMx      = (Xt0-o*Xt0(Y0,:))*(C3-C0)'+o*YM0(Y0,:);
YMg      = (Xt0-o*Xt0(Y0,:))*(C4-C0)'+o*YM0(Y0,:);
YMall    = (Xt0-o*Xt0(Y0,:))*C'+o*YM0(Y0,:);
YMnox    = (Xt0-o*Xt0(Y0,:))*(C1+C2+C4-2*C0)'+o*YM0(Y0,:);
YMnoz    = (Xt0-o*Xt0(Y0,:))*(C2+C3+C4-2*C0)'+o*YM0(Y0,:);
YMnol    = (Xt0-o*Xt0(Y0,:))*(C1+C3+C4-2*C0)'+o*YM0(Y0,:);%added
YMnog    = (Xt0-o*Xt0(Y0,:))*(C1+C2+C3-2*C0)'+o*YM0(Y0,:);%added
% the rule is to add C1-C0,C2-C0,C3-C0 or C4-C0 for each of A,L,X or G that
% is on and add C0 for each wedge that is off.
YM2zl    = (Xt0-o*Xt0(Y0,:))*(C1+C2)'+o*YM0(Y0,:);
YM2nzl   = (Xt0-o*Xt0(Y0,:))*(C3+C4)'+o*YM0(Y0,:);

YMz_gc      = (Xt0_gc-o*Xt0_gc(Y0,:))*(C1gc-C0gc)'+o*YM0_gc(Y0,:);
YMl_gc      = (Xt0_gc-o*Xt0_gc(Y0,:))*(C2gc-C0gc)'+o*YM0_gc(Y0,:);
YMx_gc      = (Xt0_gc-o*Xt0_gc(Y0,:))*(C3gc-C0gc)'+o*YM0_gc(Y0,:);
YMg_gc      = (Xt0_gc-o*Xt0_gc(Y0,:))*(C4gc-C0gc)'+o*YM0_gc(Y0,:);
YMnox_gc    = (Xt0_gc-o*Xt0_gc(Y0,:))*(C1gc+C2gc+C4gc-2*C0gc)'+o*YM0_gc(Y0,:);
YMnoz_gc    = (Xt0_gc-o*Xt0_gc(Y0,:))*(C2gc+C3gc+C4gc-2*C0gc)'+o*YM0_gc(Y0,:);
YMnol_gc    = (Xt0_gc-o*Xt0_gc(Y0,:))*(C1gc+C3gc+C4gc-2*C0gc)'+o*YM0_gc(Y0,:);%added
YMnog_gc    = (Xt0_gc-o*Xt0_gc(Y0,:))*(C1gc+C2gc+C3gc-2*C0gc)'+o*YM0_gc(Y0,:);%added


% YMn_nx     = (Xt0-o*Xt0(Y0,:))*C0nx'+o*YM0(Y0,:);
% YMz_nx      = (Xt0-o*Xt0(Y0,:))*(C1nx-C0nx)'+o*YM0(Y0,:);
% YMl_nx      = (Xt0-o*Xt0(Y0,:))*(C2nx-C0nx)'+o*YM0(Y0,:);
% YMx_nx      = (Xt0-o*Xt0(Y0,:))*(C3nx-C0nx)'+o*YM0(Y0,:);
% YMg_nx      = (Xt0-o*Xt0(Y0,:))*(C4nx-C0nx)'+o*YM0(Y0,:);
% YMall_nx    = (Xt0-o*Xt0(Y0,:))*C'+o*YM0(Y0,:);
% YMnox_nx    = (Xt0-o*Xt0(Y0,:))*(C1nx+C2nx+C4nx-2*C0nx)'+o*YM0(Y0,:);
% YMnoz_nx    = (Xt0-o*Xt0(Y0,:))*(C2nx+C3nx+C4nx-2*C0nx)'+o*YM0(Y0,:);
% YMnol_nx    = (Xt0-o*Xt0(Y0,:))*(C1nx+C3nx+C4nx-2*C0nx)'+o*YM0(Y0,:);%added
% YMnog_nx    = (Xt0-o*Xt0(Y0,:))*(C1nx+C2nx+C3nx-2*C0nx)'+o*YM0(Y0,:);%added

%consumption with only one wedge
Ymn_C=exp(YMn(:,1))-exp(YMn(:,2))-exp(lg);
Ymz_C=exp(YMz(:,1))-exp(YMz(:,2))-exp(lg);
Yml_C=exp(YMl(:,1))-exp(YMl(:,2))-exp(lg);
Ymx_C=exp(YMx(:,1))-exp(YMx(:,2))-exp(lg);
Ymg_C=exp(YMg(:,1))-exp(YMg(:,2))-exp(lgt);

%consumption with all but one wedge
Ymnoz_C=exp(YMnoz(:,1))-exp(YMnoz(:,2))-exp(lgt);
Ymnol_C=exp(YMnol(:,1))-exp(YMnol(:,2))-exp(lgt);
Ymnox_C=exp(YMnox(:,1))-exp(YMnox(:,2))-exp(lgt);
Ymnog_C=exp(YMnog(:,1))-exp(YMnog(:,2))-exp(lg);

%consumption with efficiency and labor wedge
YM2zl_C=exp(YM2zl(:,1))-exp(YM2zl(:,2))-exp(lg);
YM2nzl_C=exp(YM2nzl(:,1))-exp(YM2nzl(:,2))-exp(lgt);

             % output             %labor              %investment     %cons
data   = [t,[exp(Y(:,1)-Y(Y0,1)),exp(Y(:,3)-Y(Y0,3)),exp(Y(:,2)-Y(Y0,2)),Ct/Ct(Y0)]*100];     


%models with only one wedge
mn= [t,[exp(YMn(:,1)-YMn(Y0,1)),exp(YMn(:,3)-YMn(Y0,3)),exp(YMn(:,2)-YMn(Y0,2)),(Ymn_C/Ymn_C(Y0))]*100];
mz= [t,[exp(YMz(:,1)-YMz(Y0,1)),exp(YMz(:,3)-YMz(Y0,3)),exp(YMz(:,2)-YMz(Y0,2)),(Ymz_C/Ymz_C(Y0))]*100];
ml= [t,[exp(YMl(:,1)-YMl(Y0,1)),exp(YMl(:,3)-YMl(Y0,3)),exp(YMl(:,2)-YMl(Y0,2)),(Yml_C/Yml_C(Y0))]*100];
mx= [t,[exp(YMx(:,1)-YMx(Y0,1)),exp(YMx(:,3)-YMx(Y0,3)),exp(YMx(:,2)-YMx(Y0,2)),(Ymx_C/Ymx_C(Y0))]*100];
mg= [t,[exp(YMg(:,1)-YMg(Y0,1)),exp(YMg(:,3)-YMg(Y0,3)),exp(YMg(:,2)-YMg(Y0,2)),(Ymg_C/Ymg_C(Y0))]*100];

%models with all but one wedge
mnoz = [t,[exp(YMnoz(:,1)-YMnoz(Y0,1)),exp(YMnoz(:,3)-YMnoz(Y0,3)), ...
              exp(YMnoz(:,2)-YMnoz(Y0,2)),(Ymnoz_C/Ymnoz_C(Y0))]*100];
mnol = [t,[exp(YMnol(:,1)-YMnol(Y0,1)),exp(YMnol(:,3)-YMnol(Y0,3)), ...
              exp(YMnol(:,2)-YMnoz(Y0,2)),(Ymnol_C/Ymnol_C(Y0))]*100];
mnox = [t,[exp(YMnox(:,1)-YMnox(Y0,1)),exp(YMnox(:,3)-YMnox(Y0,3)), ...
              exp(YMnox(:,2)-YMnoz(Y0,2)),(Ymnox_C/Ymnox_C(Y0))]*100];
mnog = [t,[exp(YMnog(:,1)-YMnog(Y0,1)),exp(YMnog(:,3)-YMnog(Y0,3)), ...
              exp(YMnog(:,2)-YMnoz(Y0,2)),(Ymnog_C/Ymnog_C(Y0))]*100];

%model with efficiency and labor wedge and with no At or tault
m2zl = [t,[exp(YM2zl(:,1)-YM2zl(Y0,1)),exp(YM2zl(:,3)-YM2zl(Y0,3)),exp(YM2zl(:,2)-YM2zl(Y0,2)),(YM2zl_C/YM2zl_C(Y0))]*100];
m2nzl= [t,[exp(YM2nzl(:,1)-YM2nzl(Y0,1)),exp(YM2nzl(:,3)-YM2nzl(Y0,3)),exp(YM2nzl(:,2)-YM2nzl(Y0,2)),(YM2nzl_C/YM2nzl_C(Y0))]*100];
          
w.mle  = mle;
w.Xt0=Xt0;
w.yt = exp(lyt-lyt(Y0));
w.ht = exp(llt-llt(Y0));
w.xt = exp(lxt - lxt(Y0));
w.gt = exp(lgt-lgt(Y0));
w.ct = Ct/Ct(Y0);
w.zt   = (Zt/Zt(Y0)).^(1-theta); 
w.tault = (1-Tault)/(1-Tault(Y0));  
w.tauxt = (1+tauxt(Y0))*ones(T,1)./(1+tauxt);
w.zbar = z; w.taulbar = taul; w.tauxbar = taux; w.gbar = g;
w.Yd    = Y; %detrended 4 observables
w.Y     = data; %t + normalized to Y(0) detrended first three observables

w.mny = mn(:,2); w.mnh=mn(:,3); w.mnx = mn(:,4); w.mnc = mn(:,5);
w.mzy = mz(:,2); w.mzh=mz(:,3); w.mzx = mz(:,4); w.mzc = mz(:,5);
w.mly = ml(:,2); w.mlh=ml(:,3); w.mlx = ml(:,4); w.mlc = ml(:,5);
w.mxy = mx(:,2); w.mxh=mx(:,3); w.mxx = mx(:,4); w.mxc = mx(:,5);
w.mgy = mg(:,2); w.mgh=mg(:,3); w.mgx = mg(:,4); w.mgc = mg(:,5);
w.mnozy = mnoz(:,2); w.mnozh=mnoz(:,3); w.mnozx = mnoz(:,4); w.mnozc = mnoz(:,5);
w.mnoly = mnol(:,2); w.mnolh=mnol(:,3); w.mnolx = mnol(:,4); w.mnolc = mnol(:,5);
w.mnoxy = mnox(:,2); w.mnoxh=mnox(:,3); w.mnoxx = mnox(:,4); w.mnoxc = mnox(:,5);
w.mnogy = mnog(:,2); w.mnogh=mnog(:,3); w.mnogx = mnog(:,4); w.mnogc = mnog(:,5);


w.m2zly = m2zl(:,2); w.m2zlh = m2zl(:,3); w.m2zlx = m2zl(:,4); w.m2zlc = m2zl(:,5);
w.m2nzly = m2nzl(:,2); w.m2nzlh = m2nzl(:,3); w.m2nzlx = m2nzl(:,4); w.m2nzlc = m2nzl(:,5);

w.lgtreal = lgt;

worktemp.w = w;
% worktemp.stage = 2;
% ind = worktemp.wsize(4);wsize = worktemp.wsize(3);


% Logged and HP filtered observables and wedges
[~,worktemp.lhpw] = hpfilter(log([worktemp.w.zt worktemp.w.tault worktemp.w.tauxt worktemp.w.gt]),400*worktemp.freq); % Case: G=GC+NX
[~,worktemp.lhpo] = hpfilter(log([worktemp.w.yt worktemp.w.ht worktemp.w.xt worktemp.w.gt worktemp.w.ct]),400*worktemp.freq);
[~,worktemp.lhps] = hpfilter(log([worktemp.w.mzy worktemp.w.mly worktemp.w.mxy worktemp.w.mgy]),400*worktemp.freq);
% hours components - for sim data report
[~,worktemp.lhpsh]= hpfilter(log([worktemp.w.mzh worktemp.w.mlh worktemp.w.mxh worktemp.w.mgh]),400*worktemp.freq);
% investment components - for sim data report
[~,worktemp.lhpsx]= hpfilter(log([worktemp.w.mzx worktemp.w.mlx worktemp.w.mxx worktemp.w.mgx]),400*worktemp.freq);
% consumption components - for sim data report
[~,worktemp.lhpsc]= hpfilter(log([worktemp.w.mzc worktemp.w.mlc worktemp.w.mxc worktemp.w.mgc]),400*worktemp.freq);

% Table II A
yzrelstd     = std(worktemp.lhpw(:,1))/std(worktemp.lhpo(:,1)); % Case: G=GC+NX
ytaulrelstd  = std(worktemp.lhpw(:,2))/std(worktemp.lhpo(:,1));
ytauxrelstd  = std(worktemp.lhpw(:,3))/std(worktemp.lhpo(:,1));
ygrelstd     = std(worktemp.lhpw(:,4))/std(worktemp.lhpo(:,1));
worktemp.tableIIA1 = [yzrelstd; ytaulrelstd; ytauxrelstd; ygrelstd]; 

yzxcorr       = xcorr(worktemp.lhpw(:,1),worktemp.lhpo(:,1),4,'Coef');
ytaulxcorr    = xcorr(worktemp.lhpw(:,2),worktemp.lhpo(:,1),4,'Coef');
ytauxcorr    = xcorr(worktemp.lhpw(:,3),worktemp.lhpo(:,1),4,'Coef');
ygxcorr       = xcorr(worktemp.lhpw(:,4),worktemp.lhpo(:,1),4,'Coef');
worktemp.tableIIA2 = [yzxcorr'; ytaulxcorr'; ytauxcorr'; ygxcorr'];
 
% Table II A - observables
yyrelstd     = std(worktemp.lhpo(:,1))/std(worktemp.lhpo(:,1)); % Case: G=GC+NX
yhrelstd     = std(worktemp.lhpo(:,2))/std(worktemp.lhpo(:,1));
yxrelstd     = std(worktemp.lhpo(:,3))/std(worktemp.lhpo(:,1));
ygtrelstd    = std(worktemp.lhpo(:,4))/std(worktemp.lhpo(:,1));
worktemp.tableIIA1o = [yyrelstd; yhrelstd; yxrelstd; ygtrelstd]; 

yyxcorr       = xcorr(worktemp.lhpo(:,1),worktemp.lhpo(:,1),4,'Coef');
yhxcorr       = xcorr(worktemp.lhpo(:,2),worktemp.lhpo(:,1),4,'Coef');
yxxcorr       = xcorr(worktemp.lhpo(:,3),worktemp.lhpo(:,1),4,'Coef');
ygtxcorr      = xcorr(worktemp.lhpo(:,4),worktemp.lhpo(:,1),4,'Coef');
worktemp.tableIIA2o = [yyxcorr'; yhxcorr'; yxxcorr'; ygtxcorr']; 

% Table II B
ztaulxcorr    = xcorr(worktemp.lhpw(:,1),worktemp.lhpw(:,2),4,'Coef'); % Case: G=GC+NX
ztauxxcorr    = xcorr(worktemp.lhpw(:,1),worktemp.lhpw(:,3),4,'Coef');
zgxcorr       = xcorr(worktemp.lhpw(:,1),worktemp.lhpw(:,4),4,'Coef');
tauxtaulxcorr = xcorr(worktemp.lhpw(:,2),worktemp.lhpw(:,3),4,'Coef');
taulgxcorr    = xcorr(worktemp.lhpw(:,2),worktemp.lhpw(:,4),4,'Coef');
tauxgxcorr    = xcorr(worktemp.lhpw(:,3),worktemp.lhpw(:,4),4,'Coef');
worktemp.tableIIB = [ztaulxcorr';ztauxxcorr';zgxcorr';tauxtaulxcorr';...
    taulgxcorr';tauxgxcorr'];

% Table II B - observables
yhxcorr       = xcorr(worktemp.lhpo(:,1),worktemp.lhpo(:,2),4,'Coef'); % Case: G=GC+NX
yxxcorr       = xcorr(worktemp.lhpo(:,1),worktemp.lhpo(:,3),4,'Coef');
ygtxcorr      = xcorr(worktemp.lhpo(:,1),worktemp.lhpo(:,4),4,'Coef');
hxxcorr       = xcorr(worktemp.lhpo(:,2),worktemp.lhpo(:,3),4,'Coef');
hgtxcorr      = xcorr(worktemp.lhpo(:,2),worktemp.lhpo(:,4),4,'Coef');
xgtxcorr      = xcorr(worktemp.lhpo(:,3),worktemp.lhpo(:,4),4,'Coef');
worktemp.tableIIBo = [yhxcorr';yxxcorr';ygtxcorr';hxxcorr';...
    hgtxcorr';xgtxcorr'];

% Table III A
yzsrelstd     = std(worktemp.lhps(:,1))/std(worktemp.lhpo(:,1)); % Case: G=GC+NX
ytaulsrelstd  = std(worktemp.lhps(:,2))/std(worktemp.lhpo(:,1));
ytauxsrelstd  = std(worktemp.lhps(:,3))/std(worktemp.lhpo(:,1));
ygsrelstd     = std(worktemp.lhps(:,4))/std(worktemp.lhpo(:,1));
worktemp.tableIIIA1 = [yzsrelstd;ytaulsrelstd;ytauxsrelstd;ygsrelstd];

yzsxcorr      = xcorr(worktemp.lhps(:,1),worktemp.lhpo(:,1),4,'Coef');
ytaulsxcorr   = xcorr(worktemp.lhps(:,2),worktemp.lhpo(:,1),4,'Coef');
ytauxsxcorr   = xcorr(worktemp.lhps(:,3),worktemp.lhpo(:,1),4,'Coef');
ygsxcorr      = xcorr(worktemp.lhps(:,4),worktemp.lhpo(:,1),4,'Coef');
worktemp.tableIIIA2 = [yzsxcorr'; ytaulsxcorr'; ytauxsxcorr'; ygsxcorr'];

% Table III A - for labor instead of output
hzsrelstd     = std(worktemp.lhpsh(:,1))/std(worktemp.lhpo(:,2)); % Case: G=GC+NX
htaulsrelstd  = std(worktemp.lhpsh(:,2))/std(worktemp.lhpo(:,2));
htauxsrelstd  = std(worktemp.lhpsh(:,3))/std(worktemp.lhpo(:,2));
hgsrelstd     = std(worktemp.lhpsh(:,4))/std(worktemp.lhpo(:,2));
worktemp.tableIIIA1h = [hzsrelstd;htaulsrelstd;htauxsrelstd;hgsrelstd];

hzsxcorr      = xcorr(worktemp.lhpsh(:,1),worktemp.lhpo(:,2),4,'Coef');
htaulsxcorr   = xcorr(worktemp.lhpsh(:,2),worktemp.lhpo(:,2),4,'Coef');
htauxsxcorr   = xcorr(worktemp.lhpsh(:,3),worktemp.lhpo(:,2),4,'Coef');
hgsxcorr      = xcorr(worktemp.lhpsh(:,4),worktemp.lhpo(:,2),4,'Coef');
worktemp.tableIIIA2h = [hzsxcorr'; htaulsxcorr'; htauxsxcorr'; hgsxcorr'];

% Table III A - for investment instead of output
xzsrelstd     = std(worktemp.lhpsx(:,1))/std(worktemp.lhpo(:,3)); % Case: G=GC+NX
xtaulsrelstd  = std(worktemp.lhpsx(:,2))/std(worktemp.lhpo(:,3));
xtauxsrelstd  = std(worktemp.lhpsx(:,3))/std(worktemp.lhpo(:,3));
xgsrelstd     = std(worktemp.lhpsx(:,4))/std(worktemp.lhpo(:,3));
worktemp.tableIIIA1x = [xzsrelstd;xtaulsrelstd;xtauxsrelstd;xgsrelstd];

xzsxcorr      = xcorr(worktemp.lhpsx(:,1),worktemp.lhpo(:,3),4,'Coef');
xtaulsxcorr   = xcorr(worktemp.lhpsx(:,2),worktemp.lhpo(:,3),4,'Coef');
xtauxsxcorr   = xcorr(worktemp.lhpsx(:,3),worktemp.lhpo(:,3),4,'Coef');
xgsxcorr      = xcorr(worktemp.lhpsx(:,4),worktemp.lhpo(:,3),4,'Coef');
worktemp.tableIIIA2x = [xzsxcorr'; xtaulsxcorr'; xtauxsxcorr'; xgsxcorr'];

% Table III A - for consumption instead of output
czsrelstd     = std(worktemp.lhpsc(:,1))/std(worktemp.lhpo(:,5)); % Case: G=GC+NX
ctaulsrelstd  = std(worktemp.lhpsc(:,2))/std(worktemp.lhpo(:,5));
ctauxsrelstd  = std(worktemp.lhpsc(:,3))/std(worktemp.lhpo(:,5));
cgsrelstd     = std(worktemp.lhpsc(:,4))/std(worktemp.lhpo(:,5));
worktemp.tableIIIA1c = [czsrelstd;ctaulsrelstd;ctauxsrelstd;cgsrelstd];

czsxcorr      = xcorr(worktemp.lhpsc(:,1),worktemp.lhpo(:,5),4,'Coef');
ctaulsxcorr   = xcorr(worktemp.lhpsc(:,2),worktemp.lhpo(:,5),4,'Coef');
ctauxsxcorr   = xcorr(worktemp.lhpsc(:,3),worktemp.lhpo(:,5),4,'Coef');
cgsxcorr      = xcorr(worktemp.lhpsc(:,4),worktemp.lhpo(:,5),4,'Coef');
worktemp.tableIIIA2c = [czsxcorr'; ctaulsxcorr'; ctauxsxcorr'; cgsxcorr'];

% Table III B
ztaulswcorr    = xcorr(worktemp.lhps(:,1),worktemp.lhpw(:,2),4,'Coef');
ztauxsxcorr    = xcorr(worktemp.lhps(:,1),worktemp.lhpw(:,3),4,'Coef');
zgswcorr       = xcorr(worktemp.lhps(:,1),worktemp.lhpw(:,4),4,'Coef');
tauxtaulscorr  = xcorr(worktemp.lhps(:,2),worktemp.lhpw(:,3),4,'Coef');
taulgscorr     = xcorr(worktemp.lhps(:,2),worktemp.lhpw(:,4),4,'Coef');
tauxgcorr      = xcorr(worktemp.lhps(:,3),worktemp.lhpw(:,4),4,'Coef');
worktemp.tableIIIB = [ztaulswcorr';ztauxsxcorr';zgswcorr';tauxtaulscorr';...
    taulgscorr';tauxgcorr'];

% Table III B - labor components
ztaulswcorr    = xcorr(worktemp.lhpsh(:,1),worktemp.lhpw(:,2),4,'Coef');
ztauxsxcorr    = xcorr(worktemp.lhpsh(:,1),worktemp.lhpw(:,3),4,'Coef');
zgswcorr       = xcorr(worktemp.lhpsh(:,1),worktemp.lhpw(:,4),4,'Coef');
tauxtaulscorr  = xcorr(worktemp.lhpsh(:,2),worktemp.lhpw(:,3),4,'Coef');
taulgscorr     = xcorr(worktemp.lhpsh(:,2),worktemp.lhpw(:,4),4,'Coef');
tauxgcorr      = xcorr(worktemp.lhpsh(:,3),worktemp.lhpw(:,4),4,'Coef');
worktemp.tableIIIBh = [ztaulswcorr';ztauxsxcorr';zgswcorr';tauxtaulscorr';...
    taulgscorr';tauxgcorr'];

% Table III B - investment components
ztaulswcorr    = xcorr(worktemp.lhpsx(:,1),worktemp.lhpw(:,2),4,'Coef');
ztauxsxcorr    = xcorr(worktemp.lhpsx(:,1),worktemp.lhpw(:,3),4,'Coef');
zgswcorr       = xcorr(worktemp.lhpsx(:,1),worktemp.lhpw(:,4),4,'Coef');
tauxtaulscorr  = xcorr(worktemp.lhpsx(:,2),worktemp.lhpw(:,3),4,'Coef');
taulgscorr     = xcorr(worktemp.lhpsx(:,2),worktemp.lhpw(:,4),4,'Coef');
tauxgcorr      = xcorr(worktemp.lhpsx(:,3),worktemp.lhpw(:,4),4,'Coef');
worktemp.tableIIIBx = [ztaulswcorr';ztauxsxcorr';zgswcorr';tauxtaulscorr';...
    taulgscorr';tauxgcorr'];

ind = worktemp.bind;
wsize = size(worktemp.time(worktemp.bind:worktemp.wend),1);
% errors for 1 wedge economies - output
    mzye =   worktemp.w.Y(ind:ind+wsize-1,2)./worktemp.w.Y(ind,2) - ...
           worktemp.w.mzy(ind:ind+wsize-1)./worktemp.w.mzy(ind);
    mlye =   worktemp.w.Y(ind:ind+wsize-1,2)./worktemp.w.Y(ind,2) - ...
           worktemp.w.mly(ind:ind+wsize-1)./worktemp.w.mly(ind);
    mxye =   worktemp.w.Y(ind:ind+wsize-1,2)./worktemp.w.Y(ind,2) - ...
           worktemp.w.mxy(ind:ind+wsize-1)./worktemp.w.mxy(ind);
    mgye =   worktemp.w.Y(ind:ind+wsize-1,2)./worktemp.w.Y(ind,2) - ...
           worktemp.w.mgy(ind:ind+wsize-1)./worktemp.w.mgy(ind);    
    worktemp.w.w1yerrors = [mzye mlye mxye mgye];
    temp1 = mean(worktemp.w.w1yerrors(:,1:4).^2,1);
    worktemp.w.w1yfz = (1/temp1(1))/sum(temp1.^-1);
    worktemp.w.w1yfl = (1/temp1(2))/sum(temp1.^-1);
    worktemp.w.w1yfx = (1/temp1(3))/sum(temp1.^-1);
    worktemp.w.w1yfg = (1/temp1(4))/sum(temp1.^-1);
    

    % errors for 1 wedge economies - hours
    mzhe =   worktemp.w.Y(ind:ind+wsize-1,3)./worktemp.w.Y(ind,3) - ...
           worktemp.w.mzh(ind:ind+wsize-1)./worktemp.w.mzh(ind);
    mlhe =   worktemp.w.Y(ind:ind+wsize-1,3)./worktemp.w.Y(ind,3) - ...
           worktemp.w.mlh(ind:ind+wsize-1)./worktemp.w.mlh(ind);
    mxhe =   worktemp.w.Y(ind:ind+wsize-1,3)./worktemp.w.Y(ind,3) - ...
           worktemp.w.mxh(ind:ind+wsize-1)./worktemp.w.mxh(ind);
    mghe =   worktemp.w.Y(ind:ind+wsize-1,3)./worktemp.w.Y(ind,3) - ...
           worktemp.w.mgh(ind:ind+wsize-1)./worktemp.w.mgh(ind);    
    worktemp.w.w1herrors = [mzhe mlhe mxhe mghe];
    temp2 = mean(worktemp.w.w1herrors(:,1:4).^2,1);
    worktemp.w.w1hfz = (1/temp2(1))/sum(temp2.^-1);
    worktemp.w.w1hfl = (1/temp2(2))/sum(temp2.^-1);
    worktemp.w.w1hfx = (1/temp2(3))/sum(temp2.^-1);
    worktemp.w.w1hfg = (1/temp2(4))/sum(temp2.^-1);
    
    % errors for 1 wedge economies - investment
    mzxe =   worktemp.w.Y(ind:ind+wsize-1,4)./worktemp.w.Y(ind,4) - ...
           worktemp.w.mzx(ind:ind+wsize-1)./worktemp.w.mzx(ind);
    mlxe =   worktemp.w.Y(ind:ind+wsize-1,4)./worktemp.w.Y(ind,4) - ...
           worktemp.w.mlx(ind:ind+wsize-1)./worktemp.w.mlx(ind);
    mxxe =   worktemp.w.Y(ind:ind+wsize-1,4)./worktemp.w.Y(ind,4) - ...
           worktemp.w.mxx(ind:ind+wsize-1)./worktemp.w.mxx(ind);
    mgxe =   worktemp.w.Y(ind:ind+wsize-1,4)./worktemp.w.Y(ind,4) - ...
           worktemp.w.mgx(ind:ind+wsize-1)./worktemp.w.mgx(ind);    
    worktemp.w.w1xerrors = [mzxe mlxe mxxe mgxe];
    temp3 = mean(worktemp.w.w1xerrors(:,1:4).^2,1);
    worktemp.w.w1xfz = (1/temp3(1))/sum(temp3.^-1);
    worktemp.w.w1xfl = (1/temp3(2))/sum(temp3.^-1);
    worktemp.w.w1xfx = (1/temp3(3))/sum(temp3.^-1);
    worktemp.w.w1xfg = (1/temp3(4))/sum(temp3.^-1);
    
    % errors for 1 wedge economies - consumption
    mzce =   worktemp.w.Y(ind:ind+wsize-1,5)./worktemp.w.Y(ind,5) - ...
           worktemp.w.mzc(ind:ind+wsize-1)./worktemp.w.mzc(ind);
    mlce =   worktemp.w.Y(ind:ind+wsize-1,5)./worktemp.w.Y(ind,5) - ...
           worktemp.w.mlc(ind:ind+wsize-1)./worktemp.w.mlc(ind);
    mxce =   worktemp.w.Y(ind:ind+wsize-1,5)./worktemp.w.Y(ind,5) - ...
           worktemp.w.mxc(ind:ind+wsize-1)./worktemp.w.mxc(ind);
    mgce =   worktemp.w.Y(ind:ind+wsize-1,5)./worktemp.w.Y(ind,5) - ...
           worktemp.w.mgc(ind:ind+wsize-1)./worktemp.w.mgc(ind);    
    worktemp.w.w1cerrors = [mzce mlce mxce mgce];
    temp31 = mean(worktemp.w.w1cerrors(:,1:4).^2,1);
    worktemp.w.w1cfz = (1/temp31(1))/sum(temp31.^-1);
    worktemp.w.w1cfl = (1/temp31(2))/sum(temp31.^-1);
    worktemp.w.w1cfx = (1/temp31(3))/sum(temp31.^-1);
    worktemp.w.w1cfg = (1/temp31(4))/sum(temp31.^-1);

save('worktemp.mat','worktemp','-mat');

