% octave_fresh_run.m
% Full from-scratch BCKM MLE estimation starting from raw data.mat.
% Replicates datamine.m setup + runmleadj.m optimization without GUI calls.
%
% Run from the matlab_reference/ directory.
% Results saved to ../octave_output/fresh_run_results.mat

clear; clc;

fprintf('=== BCKM FRESH FROM-SCRATCH MLE RUN ===\n');
fprintf('Starting from raw data.mat — no previous results used.\n\n');

% ------------------------------------------------------------------ %
% 1. Load raw data (datamine.m lines 3-6)                             %
% ------------------------------------------------------------------ %
load data.mat
t    = (1980.25:0.25:2015)';
ypc  = data(:,1);
xpc  = data(:,2);
hpc  = data(:,3);
gpc  = data(:,4);
cpc  = data(:,5);
iP   = data(:,6);

nobs = size(ypc,1);
fprintf('Data loaded: T = %d quarters (1980Q1-2014Q4)\n', nobs);

% ------------------------------------------------------------------ %
% 2. Calibration (datamine.m lines 10-15)                             %
% ------------------------------------------------------------------ %
gn    = (iP(end)/iP(1))^(1/(nobs-1)) - 1;
beta  = 0.975^(1/4);
delta = 1-(1-0.05)^(1/4);
psi   = 2.5;
sigma = 1.000001;
theta = 1/3;

fprintf('\nCalibration parameters:\n');
fprintf('  gn (quarterly pop growth) = %.6f\n', gn);
fprintf('  beta (quarterly)          = %.6f\n', beta);
fprintf('  delta (quarterly)         = %.6f\n', delta);
fprintf('  psi                       = %.4f\n',  psi);
fprintf('  sigma                     = %.6f\n', sigma);
fprintf('  theta (capital share)     = %.6f\n', theta);

% ------------------------------------------------------------------ %
% 3. Detrend via maketrend (datamine.m lines 27-31)                   %
% ------------------------------------------------------------------ %
bdate    = 2008.25;
bind     = find(t == bdate);
mlestart = 1980.25;
mleend   = 2015;
iobs     = find(t == mlestart);
eobs     = find(t == mleend);
mlep     = [iobs eobs];
nps      = 50;
pb       = 0.99;
adja     = 12.88;   % adjc=2 (BGG)

[mled, Y, gz] = maketrend(t, ypc, xpc, hpc, gpc, cpc, bind, mlep);
param = [gn; gz; beta; delta; psi; sigma; theta];

fprintf('  gz (quarterly tech growth) = %.6f  (solved by maketrend/calgz)\n', gz);

% ------------------------------------------------------------------ %
% 4. Build and save worktemp (datamine.m lines 34-58)                 %
%    We save WITHOUT any mle/mlemax fields so this is a clean slate.  %
% ------------------------------------------------------------------ %
worktemp.optimnum.nps = nps;
worktemp.optimnum.pb  = pb;
worktemp.mlestart     = mlestart;
worktemp.mleend       = mleend;
worktemp.bind         = bind;
worktemp.bdate        = bdate;
worktemp.iobs         = iobs;
worktemp.eobs         = eobs;
worktemp.wend         = find(t == 2011.75);
worktemp.time         = mled(:,1);
worktemp.mled         = mled;
worktemp.Y            = Y;
worktemp.cname        = 'USA';
worktemp.freq         = 4;
worktemp.obs          = [ypc xpc hpc gpc cpc];
worktemp.params       = param;
worktemp.adjc         = 2;

save('worktemp.mat', 'worktemp', '-mat');
fprintf('\nFresh worktemp.mat saved (no prior mle results).\n');

% ------------------------------------------------------------------ %
% 5. Solve for initial Sbar (runmleadj.m line 14)                     %
% ------------------------------------------------------------------ %
fprintf('\nSolving for initial Sbar via fsolve(@initmle)...\n');
sbari = fsolve(@initmle, [0; 0.05; 0.0; log(0.2)]);
fprintf('  sbari = [%.6f, %.6f, %.6f, %.6f]\n', ...
        sbari(1), sbari(2), sbari(3), sbari(4));

% ------------------------------------------------------------------ %
% 6. Build x0b warm-start (runmleadj.m lines 48-78, adjc=2)          %
% ------------------------------------------------------------------ %
x0 = [sbari(1); sbari(2); sbari(3); sbari(4);
      0.995; 0;    0;    0;
      0;    0.995; 0;    0;
      0;    0;    0.995; 0;
      0;    0;    0;    0.995;
      0.02382866874935; -0.01028081000637;  0.00553372214162; 0;
      0.02739908940082; -0.03354291115539;  0;
      0.04951600718821;  0;
      0.10194597840286];

fprintf('\nWarm-start x0b built (adjc=2, adja=%.2f).\n', adja);
[L_x0, ~, ~, ~, ~, ~, ~, ~, ~] = mleqadj(x0, adja);
fprintf('Objective at warm-start: F = %.6f\n\n', L_x0);

% ------------------------------------------------------------------ %
% 7. Run optimization loop (runmleadj.m lines 120-141, no waitbar)    %
% ------------------------------------------------------------------ %
fprintf('Running %d restarts (pb=%.2f)...\n\n', nps, pb);

X = zeros(30, nps);
F = zeros(nps, 1);

[x1, f, ~, ~, ~] = uncmin(x0, 'mleqadj', adja);
X(:,1) = x1;  F(1) = f;
fprintf('  Restart  1: F = %.6f\n', f);

[x1, f, ~, ~, ~] = uncmin(x1, 'mleqadj', adja);
X(:,2) = x1;  F(2) = f;
fprintf('  Restart  2: F = %.6f\n', f);

x2 = x1;
for i = 3:nps
    [x2, f, ~, ~, ~] = uncmin(x2 * pb, 'mleqadj', adja);
    F(i)   = f;
    X(:,i) = x2;
    fprintf('  Restart %2d: F = %.6f\n', i, f);
end

% ------------------------------------------------------------------ %
% 8. Select best result (runmleadj.m lines 143-148)                   %
% ------------------------------------------------------------------ %
[F_best, i_best] = min(F);
if numel(i_best) > 1
    i_best = i_best(1);
end
x_best = X(:, i_best);

% ------------------------------------------------------------------ %
% 9. Extract and report parameters                                     %
% ------------------------------------------------------------------ %
Sbar_fresh = x_best(1:4);
% P stored column-major: Theta(5:8)=col1, (9:12)=col2, (13:16)=col3, (17:20)=col4
P_fresh    = reshape(x_best(5:20), 4, 4);
P0_fresh   = (eye(4) - P_fresh) * Sbar_fresh;
eigs_fresh = abs(eig(P_fresh));

% Reconstruct Q as symmetric (same as gmle.m)
Q_fresh = zeros(4,4);
Q_fresh(1,1)=x_best(21); Q_fresh(1,2)=x_best(22); Q_fresh(1,3)=x_best(23); Q_fresh(1,4)=x_best(24);
Q_fresh(2,1)=x_best(22); Q_fresh(2,2)=x_best(25); Q_fresh(2,3)=x_best(26); Q_fresh(2,4)=x_best(27);
Q_fresh(3,1)=x_best(23); Q_fresh(3,2)=x_best(26); Q_fresh(3,3)=x_best(28); Q_fresh(3,4)=x_best(29);
Q_fresh(4,1)=x_best(24); Q_fresh(4,2)=x_best(27); Q_fresh(4,3)=x_best(29); Q_fresh(4,4)=x_best(30);
V_fresh = Q_fresh * Q_fresh';

% Published values for comparison
P0_pub  = [0.013982; 0.000787; 0.012883; -0.013697];
Pdiag_pub = [0.988692; 1.001057; 0.967506; 0.994461];

fprintf('\n=== FRESH ESTIMATION RESULTS ===\n\n');
fprintf('Best restart: %d of %d\n', i_best, nps);
fprintf('Log-likelihood (-F_best): %.6f\n', -F_best);
fprintf('Published LL (worktemp):  %.6f\n', 2402.876124);
fprintf('Gap:                      %.6f\n\n', 2402.876124 - (-F_best));

fprintf('           Sbar     P_diag      P0\n');
fprintf('         fresh pub  fresh pub   fresh    pub\n');
labs = {'A ','tL','tX','g '};
Pdiag_pub2 = [0.988692; 1.001057; 0.967506; 0.994461];
Sbar_pub   = [0.133605; 0.369144; -0.046008; -1.935515];
for j = 1:4
    fprintf('  %s:  %7.4f %7.4f  %7.4f %7.4f  %8.5f %8.5f\n', ...
        labs{j}, Sbar_fresh(j), Sbar_pub(j), ...
        P_fresh(j,j), Pdiag_pub2(j), ...
        P0_fresh(j), P0_pub(j));
end

fprintf('\nMax |eig(P)| fresh: %.6f\n', max(eigs_fresh));
fprintf('Max |eig(P)| pub:   %.6f\n', 0.995133);

% ------------------------------------------------------------------ %
% 10. Save                                                             %
% ------------------------------------------------------------------ %
worktemp.mlemax.F      = F;
worktemp.mlemax.X      = X;
worktemp.mlemax.X0     = x0;
worktemp.mlemax.nps    = nps;
worktemp.mlemax.optind = i_best;
save('worktemp.mat', 'worktemp', '-mat');

outfile = '../octave_output/fresh_run_results.mat';
save(outfile, 'F', 'X', 'x_best', 'F_best', 'i_best', ...
     'Sbar_fresh', 'P_fresh', 'P0_fresh', 'Q_fresh', 'V_fresh', ...
     'param', 'gz', '-mat');
fprintf('\nResults saved to %s\n', outfile);
