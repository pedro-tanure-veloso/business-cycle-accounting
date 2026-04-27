function xpto = datamine()
t   = (1980.25:0.25:2015)';

load data.mat
ypc = data(:,1); xpc = data(:,2); hpc = data(:,3); gpc = data(:,4);
cpc = data(:,5); iP  = data(:,6);

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
worktemp.optimnum.nps = nps;
worktemp.optimnum.pb  = pb;
worktemp.mlestart     = mlestart;
worktemp.mleend       = mleend;
worktemp.bind         = bind;
worktemp.bdate        = bdate;
worktemp.iobs         = iobs;
worktemp.eobs         = eobs;
worktemp.wend         = find(t==2011.75); %event window is worktemp.bind to worktemp.wend

% data and parametrization
worktemp.time         = mled(:,1);
worktemp.mled         = mled;
worktemp.Y            = Y;
worktemp.cname        ='USA';
worktemp.freq         =4;      % number of observations per year
worktemp.obs          =[ypc xpc hpc gpc cpc];
worktemp.params       = param;
worktemp.adjc         = 2;     %1 for no, 2 for BGG, 3 for 4*BGG.   

%% save and initialize gmle.m

save('worktemp.mat','worktemp','-mat');
run('gmle.m');run('gwedges2.m');



