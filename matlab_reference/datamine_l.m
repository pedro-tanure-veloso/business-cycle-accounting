function xpto = datamine()
t   = (1980.25:0.25:2015)';

load data_l.mat
ypc = data_l(:,1); xpc = data_l(:,2); hpc = data_l(:,3); gpc = data_l(:,4);
cpc = data_l(:,5); iP  = data_l(:,6);

% choose the calibration values
nobs  = size(ypc,1);
gn    = (iP(end)/iP(1))^(1/(size(ypc,1)-1))-1; 
beta  = .975^(1/4); 
delta = 1-(1-0.05)^(1/4);
psi   = 2.5;
sigma = 1.000001;
theta = 1/3;

% choose MLE and data transformation parameters
bdate    = 2008.25;         % base date, through whih we normalize the data
bind     = find(t==bdate); % find index number of base date
mlestart = 1980.25;        % set start date for the MLE estimation
mleend   = 2015;           % set end date for the MLE estimation
iobs     = find(t==mlestart); % find index number of startmle
eobs     = find(t==mleend);   % find index number of endtmle
mlep     = [iobs eobs];    % stack initial and end index of mle sample
ssize    = eobs-iobs+1;    % mle sample size
nps      = 50;             % number of uncmin runs
pb       = 0.99;           % x0[k+1] = x0[k]*pb restart point for mle 

% detrend the variables
[mled,Y,gz]=maketrend(t,ypc,xpc,hpc,gpc,cpc,bind,mlep);
param = [gn;gz;beta;delta;psi;sigma;theta];

%% create bca file structure

% mle stuff
worktemp_l.optimnum.nps = nps;
worktemp_l.optimnum.pb  = pb;
worktemp_l.mlestart     = mlestart;
worktemp_l.mleend       = mleend;
worktemp_l.bind         = bind;
worktemp_l.bdate        = bdate;
worktemp_l.iobs         = iobs;
worktemp_l.eobs         = eobs;
worktemp_l.wend         = find(t==2011.75); %event window is worktemp.bind to worktemp.wend

% data and parametrization
worktemp_l.time         = mled(:,1);
worktemp_l.mled         = mled;
worktemp_l.Y            = Y;
worktemp_l.cname        ='USA';
worktemp_l.freq         =4;      % number of observations per year
worktemp_l.obs          =[ypc xpc hpc gpc cpc];
worktemp_l.params       = param;
worktemp_l.adjc         = 2;     %1 for no, 2 for BGG, 3 for 4*BGG.   

%% save and initialize gmle.m

save('worktemp_l.mat','worktemp_l','-mat');
run('gmle.m');run('gwedges2.m');



