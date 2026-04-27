function mle = gmle()

% Assumes data folder to be current folder.

% loads country input structure file
load('worktemp.mat');

% declares/loads observables and parametrization as global variables for
% runmle.m, mleq.m and mleseq.m

% starts MLE estimation
[Theta,Lk] = runmleadj(); load worktemp.mat

% organizes output for Sbar, P0 and P
sbar = zeros(4,1)*NaN;P=zeros(4,4)*NaN;Q=zeros(4,4)*NaN;

Sbar(1,1)=Theta(1); Sbar(2,1)=Theta(2); 
Sbar(3,1)=Theta(3); Sbar(4,1)=Theta(4);

P(1,1) = Theta(5); P(2,1) = Theta(6); P(3,1) = Theta(7); P(4,1) = Theta(8);
P(1,2) = Theta(9); P(2,2) = Theta(10);P(3,2) = Theta(11);P(4,2) = Theta(12);
P(1,3) = Theta(13);P(2,3) = Theta(14);P(3,3) = Theta(15);P(4,3) = Theta(16);
P(1,4) = Theta(17);P(2,4) = Theta(18);P(3,4) = Theta(19);P(4,4) = Theta(20);
Q(1,1) = Theta(21);Q(1,2) = Theta(22);Q(1,3) = Theta(23);Q(1,4) = Theta(24);
Q(2,1) = Q(1,2)   ;Q(2,2) = Theta(25);Q(2,3) = Theta(26);Q(2,4) = Theta(27);
Q(3,1) = Q(1,3)   ;Q(3,2) = Q(2,3)   ;Q(3,3) = Theta(28);Q(3,4) = Theta(29);
Q(4,1) = Q(1,4)   ;Q(4,2) = Q(2,4)   ;Q(4,3) = Q(3,4)   ;Q(4,4) = Theta(30);

P0      = (eye(4)-P)*Sbar;

% %organizes the standard errors for the above objects
% sesbar = zeros(4,1)*NaN;seP=zeros(4,4)*NaN;seQ=zeros(4,4)*NaN;
% 
% seP(1,1) = ses(5); seP(2,1) = ses(6); seP(3,1) = ses(7); seP(4,1) = ses(8);
% seP(1,2) = ses(9); seP(2,2) = ses(10);seP(3,2) = ses(11);seP(4,2) = ses(12);
% seP(1,3) = ses(13);seP(2,3) = ses(14);seP(3,3) = ses(15);seP(4,3) = ses(16);
% seP(1,4) = ses(17);seP(2,4) = ses(18);seP(3,4) = ses(19);seP(4,4) = ses(20);
% seQ(1,1) = ses(21);seQ(2,1) = ses(22);seQ(3,1) = ses(23);seQ(4,1) = ses(24);
% seQ(2,2) = ses(25);seQ(3,2) = ses(26);seQ(4,2) = ses(27);seQ(3,3) = ses(28);
% seQ(4,3) = ses(29);seQ(4,4) = ses(30);

% creates structure variable with all estimation outputs and inputs

mle.Theta = Theta;
mle.sbar.estimate = Sbar; %mle.sbar.se = sesbar;
mle.P.estimate = P; %mle.P.se = seP;
mle.Q.estimate = Q; %mle.Q.se = seQ;
mle.P0 = P0;
mle.Likelihood = Lk;
mle.obs = worktemp.mled(:,2:6);
mle.params = worktemp.params;
% mle.dtvars = worktemp.Y;
worktemp.mle = mle;

% saves structure variable with country name
save('worktemp.mat','worktemp','-mat');










