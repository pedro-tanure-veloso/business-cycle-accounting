% octave_full_pipeline.m
%
% Runs the complete BCKM counterfactual + f-stat pipeline for the USA
% using the published MLE solution, then exports all outputs to
% ../octave_output/ for comparison with the Python pipeline.
%
% What this script does:
%   1. Backs up worktemp.mat (restores it at the end — no permanent changes)
%   2. Ensures worktemp.mle holds the published optimal theta
%      (injects it from the hard-coded values in STEP1_OCTAVE_RESULTS.md
%       if worktemp.mle is missing or has a worse likelihood)
%   3. Runs gwedges2.m  →  wedge paths, counterfactual paths, Tables IIA/IIB/IIIA
%   4. Exports CSVs and a human-readable results file to ../octave_output/
%   5. Restores original worktemp.mat
%
% Does NOT re-run MLE estimation (that takes ~20 min; already done in
% octave_fresh_run.m and verified in STEP1_OCTAVE_RESULTS.md).
%
% Usage:
%   cd matlab_reference
%   octave --no-gui octave_full_pipeline.m 2>&1 | tee ../octave_output/full_pipeline_log.txt

fprintf('=== BCKM Full Pipeline Replication ===\n');
fprintf('Date: %s\n\n', datestr(now));

% ---------------------------------------------------------------------- %
% 0.  Directories and backup                                              %
% ---------------------------------------------------------------------- %
outdir = '../octave_output';
if ~exist(outdir, 'dir'); mkdir(outdir); end

fprintf('Backing up worktemp.mat -> %s/worktemp_pre_full_pipeline.mat\n', outdir);
copyfile('worktemp.mat', [outdir '/worktemp_pre_full_pipeline.mat']);

% ---------------------------------------------------------------------- %
% 1.  Load worktemp and validate / inject published MLE solution          %
% ---------------------------------------------------------------------- %
load worktemp.mat

% Published optimal theta from STEP1_OCTAVE_RESULTS.md / bckm_summary.txt
% F_pub = -2402.876124  (i.e. LL_pub = +2402.876124)
theta_pub = [ ...
   0.133605;   0.369144;  -0.0460083; -1.93552;    % Sbar(1..4)
   0.988692;   0.0306875; -0.00885811;-0.0407158;  % P col 1
  -0.0011749;  1.00106;   -0.0275264;  0.0174932;  % P col 2
  -0.00449176; 0.0448986;  0.967506;  -0.0425751;  % P col 3
   0.00632591; 0.00165897; 0.00156723; 0.994461;   % P col 4
   0.00768662; 0.00237971;-0.00411471; 0.000332668;% Q row 1
   0.00425474; 0.00229478; 0.0153005;              % Q row 2 (lower tri)
   0.00884692; 0.0121407;                          % Q row 3
   0.0139378];                                     % Q(4,4)

F_pub = -2402.876124;

need_inject = true;
if isfield(worktemp, 'mle') && isfield(worktemp.mle, 'Likelihood')
    if worktemp.mle.Likelihood <= F_pub + 0.1   % within 0.1 nats -> accept
        fprintf('worktemp.mle already present (L=%.6f). Using it.\n', ...
                worktemp.mle.Likelihood);
        need_inject = false;
    else
        fprintf('worktemp.mle.Likelihood=%.6f worse than published. Re-injecting.\n', ...
                worktemp.mle.Likelihood);
    end
end

if need_inject
    fprintf('Injecting published MLE solution (L=%.6f)...\n', -F_pub);

    Theta = theta_pub;
    Sbar  = Theta(1:4);
    P     = reshape(Theta(5:20), 4, 4);
    % Q is symmetric: upper-triangle storage in theta(21:30)
    Q = zeros(4,4);
    Q(1,1)=Theta(21); Q(1,2)=Theta(22); Q(1,3)=Theta(23); Q(1,4)=Theta(24);
    Q(2,1)=Theta(22); Q(2,2)=Theta(25); Q(2,3)=Theta(26); Q(2,4)=Theta(27);
    Q(3,1)=Theta(23); Q(3,2)=Theta(26); Q(3,3)=Theta(28); Q(3,4)=Theta(29);
    Q(4,1)=Theta(24); Q(4,2)=Theta(27); Q(4,3)=Theta(29); Q(4,4)=Theta(30);
    P0    = (eye(4) - P) * Sbar;

    mle.Theta      = Theta;
    mle.sbar.estimate = Sbar;
    mle.P.estimate = P;
    mle.Q.estimate = Q;
    mle.P0         = P0;
    mle.Likelihood = F_pub;             % stored as negative F (i.e. -LL)
    mle.obs        = worktemp.mled(:, 2:6);
    mle.params     = worktemp.params;

    worktemp.mle = mle;
    save('worktemp.mat', 'worktemp', '-mat');
    fprintf('Injected and saved.\n\n');
end

% ---------------------------------------------------------------------- %
% 2.  Verify mleqadj at the published theta                               %
% ---------------------------------------------------------------------- %
adja = 12.88;   % adjc=2 (BGG), matches datamine.m adjcs(2)
fprintf('Verifying mleqadj at published theta...\n');
[F_check, ~, ~, ~, ~, ~, ~, ~, ~] = mleqadj(theta_pub, adja);
fprintf('  F_check = %.6f  (published F = %.6f, gap = %.2e)\n\n', ...
        F_check, F_pub, abs(F_check - F_pub));

% ---------------------------------------------------------------------- %
% 3.  Run gwedges2.m                                                      %
%     This adds worktemp.w (wedges, counterfactual paths, Tables IIA/IIB) %
%     and saves to worktemp.mat                                            %
% ---------------------------------------------------------------------- %
fprintf('Running gwedges2.m (counterfactuals + HP tables)...\n');
run('gwedges2.m');
load worktemp.mat
fprintf('gwedges2.m complete.\n\n');

% ---------------------------------------------------------------------- %
% 4.  Print key findings for verification                                 %
% ---------------------------------------------------------------------- %
fprintf('=== KEY FINDINGS ===\n\n');
fprintf('--- Table 9: P0 ---\n');
fprintf('           fresh   BCKM Table 9\n');
P0vals = worktemp.mle.P0;
P0pub  = [0.0140; 0.0008; 0.0129; -0.0137];
labs   = {'A ', 'tL', 'tX', 'g '};
for j = 1:4
    fprintf('  %s:  %8.5f   %8.4f\n', labs{j}, P0vals(j), P0pub(j));
end

fprintf('\n--- Table 8: P diagonal ---\n');
P_est = worktemp.mle.P.estimate;
Pdiag_pub = [0.9887; 1.0011; 0.9675; 0.9945];
for j = 1:4
    fprintf('  %s:  %8.5f   %8.4f\n', labs{j}, P_est(j,j), Pdiag_pub(j));
end
eigs_P = abs(eig(P_est));
fprintf('  max|eig(P)| = %.6f\n', max(eigs_P));

fprintf('\n--- Table 11: F-statistics (1-wedge economies, output Y) ---\n');
fprintf('           ours     BCKM\n');
fprintf('  fY[A]:   %.4f   0.16\n', worktemp.w.w1yfz);
fprintf('  fY[tL]:  %.4f   0.46\n', worktemp.w.w1yfl);
fprintf('  fY[tX]:  %.4f   0.32\n', worktemp.w.w1yfx);
fprintf('  fY[g]:   %.4f   0.06\n', worktemp.w.w1yfg);

fprintf('\n--- Table 11: F-statistics (1-wedge economies, hours H) ---\n');
fprintf('  fH[A]:   %.4f   0.04\n', worktemp.w.w1hfz);
fprintf('  fH[tL]:  %.4f   0.70\n', worktemp.w.w1hfl);
fprintf('  fH[tX]:  %.4f   0.25\n', worktemp.w.w1hfx);
fprintf('  fH[g]:   %.4f   0.01\n', worktemp.w.w1hfg);

fprintf('\n--- Table 11: F-statistics (1-wedge economies, investment X) ---\n');
fprintf('  fX[A]:   %.4f   0.05\n', worktemp.w.w1xfz);
fprintf('  fX[tL]:  %.4f   0.05\n', worktemp.w.w1xfl);
fprintf('  fX[tX]:  %.4f   0.88\n', worktemp.w.w1xfx);
fprintf('  fX[g]:   %.4f   0.02\n', worktemp.w.w1xfg);

fprintf('\n--- Great Recession peak-to-trough (2008Q1=bind -> 2009Q2=bind+5) ---\n');
bind = worktemp.bind;
gr_end = bind + 5;  % 2009Q2 is bind+5 (5 quarters after 2008Q1)
if gr_end > length(worktemp.time); gr_end = length(worktemp.time); end
fprintf('  Actual  Y trough: %.4f  (normalized to 1.00 at 2008Q1)\n', worktemp.w.Y(gr_end,2));
fprintf('  CF(A)   Y trough: %.4f\n', worktemp.w.mzy(gr_end));
fprintf('  CF(tL)  Y trough: %.4f\n', worktemp.w.mly(gr_end));
fprintf('  CF(tX)  Y trough: %.4f\n', worktemp.w.mxy(gr_end));
fprintf('  CF(g)   Y trough: %.4f\n', worktemp.w.mgy(gr_end));

% ---------------------------------------------------------------------- %
% 5.  Export results to octave_output/                                    %
% ---------------------------------------------------------------------- %
fprintf('\n=== EXPORTING ===\n');
T = length(worktemp.time);

% 5a. Full worktemp struct (scipy.io.loadmat readable)
outfile_mat = [outdir '/full_pipeline_worktemp.mat'];
save(outfile_mat, 'worktemp', '-v7');
fprintf('  %s\n', outfile_mat);

% 5b. Observed + wedge paths
outfile_wedges = [outdir '/bckm_wedges.csv'];
fid = fopen(outfile_wedges, 'w');
fprintf(fid, 'time,yt,ht,xt,gt,zt,tault,tauxt\n');
for i = 1:T
    fprintf(fid, '%.4f,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g\n', ...
            worktemp.time(i), ...
            worktemp.w.yt(i), worktemp.w.ht(i), worktemp.w.xt(i), worktemp.w.gt(i), ...
            worktemp.w.zt(i), worktemp.w.tault(i), worktemp.w.tauxt(i));
end
fclose(fid);
fprintf('  %s  (%d rows)\n', outfile_wedges, T);

% 5c. Counterfactual paths - output Y (normalized to 100 at base date)
outfile_cfy = [outdir '/bckm_cf_output.csv'];
fid = fopen(outfile_cfy, 'w');
fprintf(fid, 'time,actual_y,cf_A,cf_L,cf_X,cf_G,cf_noA,cf_noL,cf_noX,cf_noG\n');
for i = 1:T
    fprintf(fid, '%.4f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f\n', ...
            worktemp.time(i), worktemp.w.Y(i,2), ...
            worktemp.w.mzy(i), worktemp.w.mly(i), worktemp.w.mxy(i), worktemp.w.mgy(i), ...
            worktemp.w.mnozy(i), worktemp.w.mnoly(i), worktemp.w.mnoxy(i), worktemp.w.mnogy(i));
end
fclose(fid);
fprintf('  %s\n', outfile_cfy);

% 5d. Counterfactual paths - hours H
outfile_cfh = [outdir '/bckm_cf_hours.csv'];
fid = fopen(outfile_cfh, 'w');
fprintf(fid, 'time,actual_h,cf_A,cf_L,cf_X,cf_G,cf_noA,cf_noL,cf_noX,cf_noG\n');
for i = 1:T
    fprintf(fid, '%.4f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f\n', ...
            worktemp.time(i), worktemp.w.Y(i,3), ...
            worktemp.w.mzh(i), worktemp.w.mlh(i), worktemp.w.mxh(i), worktemp.w.mgh(i), ...
            worktemp.w.mnozh(i), worktemp.w.mnolh(i), worktemp.w.mnoxh(i), worktemp.w.mnogh(i));
end
fclose(fid);
fprintf('  %s\n', outfile_cfh);

% 5e. Counterfactual paths - investment X
outfile_cfx = [outdir '/bckm_cf_investment.csv'];
fid = fopen(outfile_cfx, 'w');
fprintf(fid, 'time,actual_x,cf_A,cf_L,cf_X,cf_G,cf_noA,cf_noL,cf_noX,cf_noG\n');
for i = 1:T
    fprintf(fid, '%.4f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f\n', ...
            worktemp.time(i), worktemp.w.Y(i,4), ...
            worktemp.w.mzx(i), worktemp.w.mlx(i), worktemp.w.mxx(i), worktemp.w.mgx(i), ...
            worktemp.w.mnozx(i), worktemp.w.mnolx(i), worktemp.w.mnoxx(i), worktemp.w.mnogx(i));
end
fclose(fid);
fprintf('  %s\n', outfile_cfx);

% 5f. F-statistics table (all rows)
outfile_fstats = [outdir '/bckm_fstats.csv'];
fid = fopen(outfile_fstats, 'w');
fprintf(fid, 'variable,fA,fL,fX,fG\n');
fprintf(fid, 'y,%.8f,%.8f,%.8f,%.8f\n', ...
        worktemp.w.w1yfz, worktemp.w.w1yfl, worktemp.w.w1yfx, worktemp.w.w1yfg);
fprintf(fid, 'h,%.8f,%.8f,%.8f,%.8f\n', ...
        worktemp.w.w1hfz, worktemp.w.w1hfl, worktemp.w.w1hfx, worktemp.w.w1hfg);
fprintf(fid, 'x,%.8f,%.8f,%.8f,%.8f\n', ...
        worktemp.w.w1xfz, worktemp.w.w1xfl, worktemp.w.w1xfx, worktemp.w.w1xfg);
fprintf(fid, 'c,%.8f,%.8f,%.8f,%.8f\n', ...
        worktemp.w.w1cfz, worktemp.w.w1cfl, worktemp.w.w1cfx, worktemp.w.w1cfg);
fclose(fid);
fprintf('  %s\n', outfile_fstats);

% 5g. HP-filtered cyclical tables (Table IIA1, IIA2, IIB) as flat CSV
outfile_tables = [outdir '/bckm_hp_tables.csv'];
fid = fopen(outfile_tables, 'w');
fprintf(fid, 'table,row,val\n');
% IIA1: relative std of wedges
rnames = {'A','tL','tX','g'};
for j = 1:4
    fprintf(fid, 'IIA1,%s,%.8f\n', rnames{j}, worktemp.tableIIA1(j));
end
% IIA1o: relative std of observables
onames = {'y','h','x','g_obs'};
for j = 1:4
    fprintf(fid, 'IIA1o,%s,%.8f\n', onames{j}, worktemp.tableIIA1o(j));
end
% IIA2: xcorr(wedge_i, y) at lags -4..+4  (9 values per row)
for j = 1:4
    for lag = 1:9
        fprintf(fid, 'IIA2_%s,lag%d,%.8f\n', rnames{j}, lag-5, worktemp.tableIIA2(j,lag));
    end
end
% IIB: wedge-wedge xcorr (6 pairs x 9 lags)
pair_names = {'z_tL','z_tX','z_g','tL_tX','tL_g','tX_g'};
for j = 1:6
    for lag = 1:9
        fprintf(fid, 'IIB_%s,lag%d,%.8f\n', pair_names{j}, lag-5, worktemp.tableIIB(j,lag));
    end
end
fclose(fid);
fprintf('  %s\n', outfile_tables);

% 5h. Detrended observables (BCKM's mled matrix)
outfile_mled = [outdir '/bckm_mled_full.csv'];
mled = worktemp.mled;
fid  = fopen(outfile_mled, 'w');
fprintf(fid, 'time,ypc_dt,xpc_dt,hpc_dt,gpc_dt,cpc_dt\n');
for i = 1:size(mled,1)
    fprintf(fid, '%.4f,%.10g,%.10g,%.10g,%.10g,%.10g\n', ...
            mled(i,1), mled(i,2), mled(i,3), mled(i,4), mled(i,5), mled(i,6));
end
fclose(fid);
fprintf('  %s\n', outfile_mled);

% ---------------------------------------------------------------------- %
% 6.  Write human-readable results summary                                %
% ---------------------------------------------------------------------- %
outfile_txt = [outdir '/full_pipeline_results.txt'];
fid = fopen(outfile_txt, 'w');
fprintf(fid, '=== BCKM Full Pipeline Results ===\n');
fprintf(fid, 'Run date: %s\n', datestr(now));
fprintf(fid, 'MLE solution: published optimal (L=%.6f)\n\n', -F_pub);

fprintf(fid, '--- Table 8: VAR Transition Matrix P ---\n');
for r = 1:4
    fprintf(fid, '  ');
    for c = 1:4
        fprintf(fid, '%9.5f', P_est(r,c));
    end
    fprintf(fid, '\n');
end
fprintf(fid, '  max|eig(P)| = %.6f\n\n', max(eigs_P));

fprintf(fid, '--- Table 9: P0 ---\n');
P0v = worktemp.mle.P0;
for j = 1:4
    fprintf(fid, '  %s: %9.5f  (BCKM: %6.4f)\n', labs{j}, P0v(j), P0pub(j));
end
fprintf(fid, '\n');

fprintf(fid, '--- Table 10: V = QQ'' ---\n');
Q_est = worktemp.mle.Q.estimate;
V_est = Q_est * Q_est';
for r = 1:4
    fprintf(fid, '  ');
    for c = 1:4
        fprintf(fid, '%12.5e', V_est(r,c));
    end
    fprintf(fid, '\n');
end
fprintf(fid, '\n');

fprintf(fid, '--- Table 11: F-statistics ---\n');
fprintf(fid, '           fA       fL       fX       fG\n');
fprintf(fid, '  y:   %.4f   %.4f   %.4f   %.4f\n', ...
        worktemp.w.w1yfz, worktemp.w.w1yfl, worktemp.w.w1yfx, worktemp.w.w1yfg);
fprintf(fid, '  h:   %.4f   %.4f   %.4f   %.4f\n', ...
        worktemp.w.w1hfz, worktemp.w.w1hfl, worktemp.w.w1hfx, worktemp.w.w1hfg);
fprintf(fid, '  x:   %.4f   %.4f   %.4f   %.4f\n', ...
        worktemp.w.w1xfz, worktemp.w.w1xfl, worktemp.w.w1xfx, worktemp.w.w1xfg);
fprintf(fid, '  c:   %.4f   %.4f   %.4f   %.4f\n', ...
        worktemp.w.w1cfz, worktemp.w.w1cfl, worktemp.w.w1cfx, worktemp.w.w1cfg);
fprintf(fid, '\n  BCKM Table 11 targets:\n');
fprintf(fid, '  y:   0.16     0.46     0.32     0.06\n\n');

fprintf(fid, '--- Table IIA1: Relative std of wedges (HP-filtered) ---\n');
for j = 1:4
    fprintf(fid, '  std(%s)/std(y) = %.4f\n', rnames{j}, worktemp.tableIIA1(j));
end
fprintf(fid, '\n');

fprintf(fid, '--- Table IIA1o: Relative std of observables (HP-filtered) ---\n');
for j = 1:4
    fprintf(fid, '  std(%s)/std(y) = %.4f\n', onames{j}, worktemp.tableIIA1o(j));
end
fprintf(fid, '\n');

fclose(fid);
fprintf('  %s\n', outfile_txt);

% ---------------------------------------------------------------------- %
% 7.  Restore original worktemp.mat                                       %
% ---------------------------------------------------------------------- %
fprintf('\nRestoring worktemp.mat from backup...\n');
copyfile([outdir '/worktemp_pre_full_pipeline.mat'], 'worktemp.mat');
fprintf('Restored.\n\n');

fprintf('=== Full pipeline replication complete ===\n');
fprintf('All outputs in %s/\n', outdir);
