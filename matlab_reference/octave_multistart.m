% octave_multistart.m
% Independent multi-start probe of the BCKM MLE objective.
% Runs N_STARTS random perturbations of the optimal theta, re-optimizes,
% and records converged F values and P diagonals to assess multimodality.
%
% Must be run from the matlab_reference/ directory.
% Results saved to ../octave_output/multistart_results.mat

clear; clc;

load worktemp.mat

adja    = 12.88;                    % adjc=2 (BGG adjustment costs)
x_opt   = worktemp.mle.Theta(:);   % 30x1 optimal published theta
L_opt   = worktemp.mle.Likelihood; % published optimal objective value

fprintf('=== BCKM INDEPENDENT MULTI-START PROBE ===\n');
fprintf('Optimal (published) F = %.6f\n\n', L_opt);

% Verify objective at published optimal point
[L_check, ~, ~, ~, ~, ~, ~, ~, ~] = mleqadj(x_opt, adja);
fprintf('Verified at published theta: F = %.6f (diff = %.2e)\n\n', ...
        L_check, abs(L_check - L_opt));

N_STARTS = 10;
rng(20260428);   % reproducible seed

F_vals    = zeros(N_STARTS, 1);
X_conv    = zeros(30, N_STARTS);
Pdiag_all = zeros(4, N_STARTS);

for i = 1:N_STARTS
    % Perturb Sbar (theta 1-4): ±10% additive noise scaled to Sbar magnitude
    % Perturb P off-diagonal (theta 6-9,10,13-20): ±5% relative
    % Perturb P diagonal (theta 5,10,15,20): ±3% relative (near stationarity)
    % Perturb Q cholesky (theta 21-30): ±20% relative
    x_pert = x_opt;
    x_pert(1:4)   = x_opt(1:4)  + 0.10 * abs(x_opt(1:4))  .* (2*rand(4,1)-1);
    x_pert(5:20)  = x_opt(5:20) .* (1 + 0.05*(2*rand(16,1)-1));
    x_pert(21:30) = x_opt(21:30).* (1 + 0.20*(2*rand(10,1)-1));

    fprintf('Start %2d: optimizing...', i);
    [x_conv, f_conv] = uncmin(x_pert, 'mleqadj', adja);

    % Extract P diagonal from converged theta
    Pdiag = [x_conv(5); x_conv(10); x_conv(15); x_conv(20)];

    F_vals(i)      = f_conv;
    X_conv(:, i)   = x_conv;
    Pdiag_all(:,i) = Pdiag;

    fprintf(' F = %.6f   P_diag = [%.4f %.4f %.4f %.4f]\n', ...
            f_conv, Pdiag(1), Pdiag(2), Pdiag(3), Pdiag(4));
end

fprintf('\n=== SUMMARY ===\n');
fprintf('Published F   : %.6f\n', L_opt);
fprintf('Best F found  : %.6f\n', min(F_vals));
fprintf('Worst F found : %.6f\n', max(F_vals));
fprintf('Spread        : %.6f\n', max(F_vals) - min(F_vals));
fprintf('\nAll F values:\n');
disp(sort(F_vals));

% Save
outfile = '../octave_output/multistart_results.mat';
save(outfile, 'F_vals', 'X_conv', 'Pdiag_all', 'x_opt', 'L_opt', 'L_check', '-mat');
fprintf('Saved results to %s\n', outfile);
