% dump_worktemp_for_diff.m
%
% Dumps fields of worktemp.mat to CSV / .mat in octave_output/ so the
% Python pipeline on the other machine can diff its observables and
% calibration against BCKM's.
%
% Usage (run from the BCA repo root):
%   octave --no-gui matlab_reference/dump_worktemp_for_diff.m
%
% Produces in octave_output/:
%   bckm_obs_raw.csv     -- T x 6: time, ypc, xpc, hpc, gpc, cpc (raw)
%   bckm_mled.csv        -- detrended observables (whatever maketrend returns)
%   bckm_dump.mat        -- the worktemp struct, v7 format (scipy.io readable)
%   bckm_summary.txt     -- human-readable means/ratios + MLE solution if present

candidates = {'worktemp.mat', ...
              'matlab_reference/worktemp.mat', ...
              'octave_output/worktemp.mat'};
found = '';
for k = 1:numel(candidates)
    if exist(candidates{k}, 'file')
        found = candidates{k};
        break
    end
end
if isempty(found)
    error('worktemp.mat not found in . , matlab_reference/, or octave_output/');
end
fprintf('Loading %s ...\n', found);
load(found);

if ~exist('octave_output', 'dir'); mkdir('octave_output'); end

obs  = worktemp.obs;
time = worktemp.time;
T    = size(obs, 1);

% --- Raw observables CSV: time, ypc, xpc, hpc, gpc, cpc ---
fid = fopen('octave_output/bckm_obs_raw.csv', 'w');
fprintf(fid, 'time,ypc,xpc,hpc,gpc,cpc\n');
for i = 1:T
    fprintf(fid, '%.4f,%.10g,%.10g,%.10g,%.10g,%.10g\n', ...
            time(i), obs(i,1), obs(i,2), obs(i,3), obs(i,4), obs(i,5));
end
fclose(fid);

% --- Detrended observables CSV ---
mled = worktemp.mled;
[Tm, Km] = size(mled);
fid = fopen('octave_output/bckm_mled.csv', 'w');
for j = 1:Km
    if j > 1; fprintf(fid, ','); end
    fprintf(fid, 'col%d', j);
end
fprintf(fid, '\n');
for i = 1:Tm
    for j = 1:Km
        if j > 1; fprintf(fid, ','); end
        fprintf(fid, '%.10g', mled(i, j));
    end
    fprintf(fid, '\n');
end
fclose(fid);

% --- Save .mat for scipy.io.loadmat (binary precision, all fields) ---
save('octave_output/bckm_dump.mat', 'worktemp', '-v7');

% --- Human-readable summary ---
fid = fopen('octave_output/bckm_summary.txt', 'w');
fprintf(fid, '=== BCKM data summary ===\n');
fprintf(fid, 'T            = %d quarters\n', T);
fprintf(fid, 'time range   = %g to %g\n', time(1), time(end));
fprintf(fid, '\n--- Raw obs means ---\n');
fprintf(fid, 'mean(ypc)    = %.10g\n', mean(obs(:, 1)));
fprintf(fid, 'mean(xpc)    = %.10g\n', mean(obs(:, 2)));
fprintf(fid, 'mean(hpc)    = %.10g\n', mean(obs(:, 3)));
fprintf(fid, 'mean(gpc)    = %.10g\n', mean(obs(:, 4)));
fprintf(fid, 'mean(cpc)    = %.10g\n', mean(obs(:, 5)));
fprintf(fid, '\n--- Ratios (BCKM convention) ---\n');
fprintf(fid, 'mean(g)/mean(y) = %.6f\n', mean(obs(:, 4)) / mean(obs(:, 1)));
fprintf(fid, 'mean(x)/mean(y) = %.6f\n', mean(obs(:, 2)) / mean(obs(:, 1)));
fprintf(fid, 'mean(c)/mean(y) = %.6f\n', mean(obs(:, 5)) / mean(obs(:, 1)));
fprintf(fid, 'mean(h)         = %.6f\n', mean(obs(:, 3)));

p = worktemp.params;
fprintf(fid, '\n--- Calibration params (worktemp.params) ---\n');
fprintf(fid, 'gn    = %.10g\n', p(1));
fprintf(fid, 'gz    = %.10g\n', p(2));
fprintf(fid, 'beta  = %.10g\n', p(3));
fprintf(fid, 'delta = %.10g\n', p(4));
fprintf(fid, 'psi   = %.10g\n', p(5));
fprintf(fid, 'sigma = %.10g\n', p(6));
fprintf(fid, 'theta = %.10g\n', p(7));
fprintf(fid, 'adjc  = %d\n', worktemp.adjc);

% Adjustment cost details if computed by gmle.m
if isfield(worktemp, 'adja')
    fprintf(fid, 'adja  = %.10g\n', worktemp.adja);
end
if isfield(worktemp, 'adjb')
    fprintf(fid, 'adjb  = %.10g\n', worktemp.adjb);
end

% MLE solution (added by gmle.m). BCKM nests it under worktemp.mle.{...}.
% Walk one level deep generically so we don't depend on field name casing
% or layout. Anything nested deeper than one level (e.g. mle.Sbar.estimate)
% is unwrapped if it has an `estimate` subfield, otherwise its struct
% summary is printed.
fprintf(fid, '\n--- MLE solution (if present) ---\n');

if isfield(worktemp, 'mle') && isstruct(worktemp.mle)
    mle = worktemp.mle;
    mle_fields = fieldnames(mle);
    for k = 1:numel(mle_fields)
        fname = mle_fields{k};
        val = mle.(fname);

        if isstruct(val) && isfield(val, 'estimate')
            % BCKM convention: mle.<param>.estimate is the converged value
            est = val.estimate;
            if isscalar(est)
                fprintf(fid, 'mle.%-12s.estimate = %.10g\n', fname, est);
            else
                fprintf(fid, 'mle.%-12s.estimate (size %dx%d) = %s\n', ...
                        fname, size(est, 1), size(est, 2), mat2str(est, 6));
            end
        elseif isnumeric(val) && isscalar(val)
            fprintf(fid, 'mle.%-12s = %.10g\n', fname, val);
        elseif isnumeric(val)
            fprintf(fid, 'mle.%-12s (size %dx%d) = %s\n', ...
                    fname, size(val, 1), size(val, 2), mat2str(val, 6));
        elseif ischar(val)
            fprintf(fid, 'mle.%-12s = "%s"\n', fname, val);
        elseif isstruct(val)
            sub = fieldnames(val);
            fprintf(fid, 'mle.%-12s (struct with fields: %s)\n', ...
                    fname, strjoin(sub, ', '));
        else
            fprintf(fid, 'mle.%-12s (skipped, type=%s)\n', fname, class(val));
        end
    end
else
    % Fall back to top-level field names just in case
    solution_fields = {'F', 'LL', 'Likelihood', 'theta', 'Sbar', 'P0', 'P', 'Q', 'V'};
    found_any = false;
    for k = 1:numel(solution_fields)
        f = solution_fields{k};
        if isfield(worktemp, f)
            found_any = true;
            val = worktemp.(f);
            if isscalar(val)
                fprintf(fid, '%-12s = %.10g\n', f, val);
            else
                fprintf(fid, '%-12s = %s (size %dx%d)\n', f, ...
                        mat2str(val, 6), size(val, 1), size(val, 2));
            end
        end
    end
    if ~found_any
        fprintf(fid, '(no mle struct and no top-level F/LL/Sbar/P/Q/V — ');
        fprintf(fid, 'top-level fields: %s)\n', strjoin(fieldnames(worktemp), ', '));
    end
end

fclose(fid);

% --- Echo summary to stdout too ---
fprintf('\n=== Files written ===\n');
fprintf('  octave_output/bckm_obs_raw.csv  (%d rows)\n', T);
fprintf('  octave_output/bckm_mled.csv     (%dx%d)\n', Tm, Km);
fprintf('  octave_output/bckm_dump.mat\n');
fprintf('  octave_output/bckm_summary.txt\n');

fprintf('\n=== BCKM data summary ===\n');
type octave_output/bckm_summary.txt;
