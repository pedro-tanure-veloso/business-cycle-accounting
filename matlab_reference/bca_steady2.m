%BCA_Steady  Start with initial levels for gdp, consumption, investment, 
%            and hours of work and use these observations to set 
%            parameters for a business cycle model with variations
%            an efficiency wedge, a labor wedge, an investment wedge
%            and a government consumption wedge. This is a modified version
%            of the codes for the 1-sector model in:
%
%               Unmeasured Investment and the Puzzling 1990s US Boom
%               by Ellen McGrattan and Ed Prescott
%
%            See details in the online appendix (Minneapolis Fed Staff 
%            Report 395).

%            Ellen McGrattan, 11-1-05
%            Revised, ERM, 3-8-16

%---------------------------------------------------------------------

load worktemp.mat;

% 1. Observed levels
GDPs    = mean(exp(worktemp.Y(:,1)));
Xs      = mean(exp(worktemp.Y(:,2)));
Hs      = mean(exp(worktemp.Y(:,3)));
gwedges = mean(exp(worktemp.Y(:,4)));
Cs      = mean(exp(worktemp.Y(:,5)));
% GDPs       = 1;
% Cs         = .7626*GDPs;
% Xs         = .2377*GDPs;
% Hs         = .2751;
% Ks         = 3.9112*GDPs;

% 2. Fixed parameters
%
beta       = worktemp.params(3); 
eta        = worktemp.params(1); %pop growth rate
gamma      = worktemp.params(2); %tech growth rate
sigma      = 1.000001;
lwedges    = 1;
xwedges    = 1;
bhat       = beta*(1+gamma)^(-sigma);
grate      = (1+eta)*(1+gamma);

% 3. Implied parameters
%

ys         = GDPs;
xs         = Xs;
cs         = Cs;
% gwedges    = GDPs-Cs-Xs;
hs         = Hs;
delta      = 1-(1-0.05)^(1/4);%xs/ks+ 1-grate;
ks         = xs/(delta+grate-1);
theta      = 1/3;%(1-bhat*(1-delta))/(bhat*xwedges) * ks/ys;
thet1      = 1-theta;
psi        = 2.5;%lwedges*(1-theta)*(1-hs)*ys/(cs*hs);
ewedges    = ys/(ks^theta*hs^thet1);
adja       = 0.25/(sum(worktemp.params([1:2 4])));


% 4. Write out results
%
fid        = fopen('bca_params2.m','w');
fprintf(fid,'adja       = %11.8f;\n',adja);
fprintf(fid,'beta       = %11.8f;\n',beta);
fprintf(fid,'delta      = %11.8f;\n',delta);
fprintf(fid,'eta        = %11.8f;\n',eta);
fprintf(fid,'gamma      = %11.8f;\n',gamma);
fprintf(fid,'psi        = %11.8f;\n',psi);
fprintf(fid,'sigma      = %11.8f;\n',sigma);
fprintf(fid,'theta      = %11.8f;\n',theta);
fprintf(fid,'ewedges    = %11.8f;\n',ewedges);
fprintf(fid,'lwedges    = %11.8f;\n',lwedges);
fprintf(fid,'xwedges    = %11.8f;\n',xwedges);
fprintf(fid,'gwedges    = %11.8f;\n',gwedges);
fprintf(fid,'ks         = %11.8f;\n',ks);
fprintf(fid,'hs         = %11.8f;\n',hs);
fprintf(fid,'ys         = %11.8f;\n',ys);
fprintf(fid,'gdps       = %11.8f;\n',ys);
fclose(fid);
