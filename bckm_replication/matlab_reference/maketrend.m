function [mled,Y,gzt]=maketrend(t,ypc,xpc,hpc,gpc,cpc,bind,mlep)
T = size(ypc,1);

%% calibrates norm to make mean detrended log(ypc) = 0 over MLE sample.
%  assumes that base date is inside the MLE sample.
global yy by mles mlee
yy   = ypc;
by   = bind;    % base date
mles = mlep(1); % starting obs index of mle sample
mlee = mlep(2); % ending obs index of mle sample
gzt = fsolve(@calgz,0);

%%
cpci = ypc-xpc-gpc; % implied consumption as opposed to real consumption cpc
mled  = [t,ypc/ypc(by)*(1+gzt)^by,xpc/ypc(by)*(1+gzt)^by,...
    hpc,gpc/ypc(by)*(1+gzt)^by,cpc/ypc(by)*(1+gzt)^by,cpci/ypc(by)*(1+gzt)^by];

Y        = log(mled(:,2:7))-log([(1+gzt).^[0:T-1]',(1+gzt).^[0:T-1]', ...
                         ones(T,1),(1+gzt).^[0:T-1]',(1+gzt).^[0:T-1]',...
                         (1+gzt).^[0:T-1]']);
