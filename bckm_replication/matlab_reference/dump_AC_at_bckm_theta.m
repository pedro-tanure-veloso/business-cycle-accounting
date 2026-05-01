% dump_AC_at_bckm_theta.m
%
% Dumps A (6x6) and C (4x6) state-space matrices built by mleqadj.m
% at the published BCKM theta stored in worktemp.mle.Theta.
% Also dumps all supporting quantities for element-wise Python diff.
%
% Run from repo root:
%   octave --no-gui --eval "addpath('matlab_reference'); dump_AC_at_bckm_theta"
%
% Outputs land in octave_output/

outdir = 'octave_output';

load worktemp.mat

% -----------------------------------------------------------------------
% 1. Calibration parameters (same as mleqadj.m lines 20-34)
% -----------------------------------------------------------------------
param = worktemp.params;
gn    = param(1);
gz    = param(2);
beta  = param(3);
delta = param(4);
psi   = param(5);
sigma = param(6);
theta = param(7);

adjcs = [0, 12.88, 4*12.88];
adja  = adjcs(worktemp.adjc);
param = [gn; gz; beta; delta; psi; sigma; theta; adja];

fprintf('=== dump_AC_at_bckm_theta ===\n');
fprintf('adja  = %.4f  (adjc index = %d)\n', adja, worktemp.adjc);
fprintf('gn    = %.8f\n', gn);
fprintf('gz    = %.8f\n', gz);
fprintf('beta  = %.8f\n', beta);
fprintf('delta = %.8f\n', delta);
fprintf('psi   = %.4f\n', psi);
fprintf('sigma = %.4f\n', sigma);
fprintf('theta = %.8f\n', theta);

% -----------------------------------------------------------------------
% 2. Extract Theta from worktemp.mle (30 free parameters)
% -----------------------------------------------------------------------
Theta = worktemp.mle.Theta(:);
fprintf('\nTheta (30 elements):\n');
fprintf('  %.16e\n', Theta);

% -----------------------------------------------------------------------
% 3. Assign Theta to Sbar, P, Q  (same as mleqadj.m lines 28-96)
% -----------------------------------------------------------------------
P_mat  = 0.995*eye(4);
Sbar   = [log(1); 0.05; 0.0; log(0.07)];
P0_vec = (eye(4)-P_mat)*Sbar;
Q_mat  = zeros(4);
D      = zeros(4);
R      = zeros(4);

Sbar(1)     = Theta(1);
Sbar(2)     = Theta(2);
Sbar(3)     = Theta(3);
Sbar(4)     = Theta(4);
P_mat(1,1)  = Theta(5);
P_mat(2,1)  = Theta(6);
P_mat(3,1)  = Theta(7);
P_mat(4,1)  = Theta(8);
P_mat(1,2)  = Theta(9);
P_mat(2,2)  = Theta(10);
P_mat(3,2)  = Theta(11);
P_mat(4,2)  = Theta(12);
P_mat(1,3)  = Theta(13);
P_mat(2,3)  = Theta(14);
P_mat(3,3)  = Theta(15);
P_mat(4,3)  = Theta(16);
P_mat(1,4)  = Theta(17);
P_mat(2,4)  = Theta(18);
P_mat(3,4)  = Theta(19);
P_mat(4,4)  = Theta(20);
Q_mat(1,1)  = Theta(21);
Q_mat(2,1)  = Theta(22);
Q_mat(3,1)  = Theta(23);
Q_mat(4,1)  = Theta(24);
Q_mat(2,2)  = Theta(25);
Q_mat(3,2)  = Theta(26);
Q_mat(4,2)  = Theta(27);
Q_mat(3,3)  = Theta(28);
Q_mat(4,3)  = Theta(29);
Q_mat(4,4)  = Theta(30);

P0_vec = (eye(4)-P_mat)*Sbar;

fprintf('\nSbar = [%.8f, %.8f, %.8f, %.8f]\n', Sbar(1), Sbar(2), Sbar(3), Sbar(4));
fprintf('max|eig(P)| = %.8f\n', max(abs(eig(P_mat))));

% -----------------------------------------------------------------------
% 4. Steady state  (same as mleqadj.m lines 146-161)
% -----------------------------------------------------------------------
tem    = (eye(4)-P_mat)\P0_vec;
zs     = exp(tem(1));
tauls  = tem(2);
tauxs  = tem(3);
gs     = exp(tem(4));
beth   = beta*(1+gz)^(-sigma);
kls    = ((1+tauxs)*(1-beth*(1-delta))/(beth*theta))^(1/(theta-1))*zs;
% A_ss (temporary, renamed to avoid clash with state-space A below)
A_ss   = (zs/kls)^(1-theta)-(1+gz)*(1+gn)+1-delta;
B_ss   = (1-tauls)*(1-theta)*kls^theta*zs^(1-theta)/psi;
ks     = (B_ss+gs)/(A_ss+B_ss/kls);
cs     = A_ss*ks-gs;
ls     = ks/kls;
ys     = ks^theta*(zs*ls)^(1-theta);
xs     = ys-cs-gs;
X0     = [log(ks); log(zs); tauls; tauxs; log(gs); 1];
Y0     = [log(ys); log(xs); log(ls); log(gs)];

fprintf('\nSteady-state values:\n');
fprintf('  ks    = %.10f\n', ks);
fprintf('  ys    = %.10f\n', ys);
fprintf('  ls    = %.10f\n', ls);
fprintf('  xs    = %.10f\n', xs);
fprintf('  gs    = %.10f\n', gs);
fprintf('  cs    = %.10f\n', cs);
fprintf('  zs    = %.10f  (=exp(Sbar(1)) at SS: %.10f)\n', zs, exp(Sbar(1)));
fprintf('  tauls = %.10f  (Sbar(2) = %.10f)\n', tauls, Sbar(2));
fprintf('  tauxs = %.10f  (Sbar(3) = %.10f)\n', tauxs, Sbar(3));
fprintf('  kls   = %.10f\n', kls);
fprintf('  beth  = %.10f\n', beth);

% -----------------------------------------------------------------------
% 5. Numerical residuals -> capital LOM coefficients
%    (same as mleqadj.m lines 167-194)
% -----------------------------------------------------------------------
Z   = [log(ks); log(ks); log(ks); log(zs); log(zs); tauls; tauls;
       tauxs;   tauxs;   log(gs); log(gs)];
del = max(abs(Z)*1e-5, 1e-8);
dR  = zeros(11,1);
for i = 1:11
  Zp    = Z;  Zm = Z;
  Zp(i) = Z(i)+del(i);
  Zm(i) = Z(i)-del(i);
  dR(i) = (res_adjust(Zp,param)-res_adjust(Zm,param))/(2*del(i));
end

a0     = dR(1);
a1     = dR(2);
a2     = dR(3);
b0     = dR(4:2:11)';
b1     = dR(5:2:11)';
tem    = roots([a0, a1, a2]);
gammak = tem(abs(tem)<1);

fprintf('\nQuadratic roots for gamma_k: %.8f  %.8f\n', real(tem(1)), real(tem(2)));
fprintf('Selected gammak = %.10f\n', gammak);

gamma  = -((a0*gammak+a1)*eye(4)+a0*P_mat')\(b0*P_mat+b1)';
gamma0 = (1-gammak)*log(ks)-gamma'*[log(zs); tauls; tauxs; log(gs)];
Gamma  = [gammak; gamma; gamma0];

fprintf('gamma  = [%.10f, %.10f, %.10f, %.10f]\n', gamma(1), gamma(2), gamma(3), gamma(4));
fprintf('gamma0 = %.10f\n', gamma0);

% -----------------------------------------------------------------------
% 6. Analytic partial derivatives for C  (mleqadj.m lines 202-216)
% -----------------------------------------------------------------------
philh  = -(psi*ys*(1-theta)+(1-theta)*(1-tauls)*ys*(1-ls)/ls*theta + ...
           (1-theta)*(1-tauls)*ys);
philk  = (psi*ys*theta+psi*(1-delta)*ks - ...
          (1-theta)*(1-tauls)*ys*(1-ls)/ls*theta)/philh;
philz  = (psi*ys*(1-theta)-(1-theta)^2*(1-tauls)*ys*(1-ls)/ls)/philh;
phill  = ((1-theta)*(1-tauls)*ys*(1-ls)/ls*(1/(1-tauls)))/philh;
philg  = (-psi*gs)/philh;
philkp = (-psi*(1+gz)*(1+gn)*ks)/philh;
phiyk  = theta+(1-theta)*philk;
phiyz  = (1-theta)*(1+philz);
phiyl  = (1-theta)*phill;
phiyg  = (1-theta)*philg;
phiykp = (1-theta)*philkp;
phixk  = -ks/xs*(1-delta);
phixkp = ks/xs*(1+gz)*(1+gn);

fprintf('\nPartial derivatives:\n');
fprintf('  philh  = %.10f\n', philh);
fprintf('  philk  = %.10f\n', philk);
fprintf('  philz  = %.10f\n', philz);
fprintf('  phill  = %.10f\n', phill);
fprintf('  philg  = %.10f\n', philg);
fprintf('  philkp = %.10f\n', philkp);
fprintf('  phiyk  = %.10f\n', phiyk);
fprintf('  phiyz  = %.10f\n', phiyz);
fprintf('  phiyl  = %.10f\n', phiyl);
fprintf('  phiyg  = %.10f\n', phiyg);
fprintf('  phiykp = %.10f\n', phiykp);
fprintf('  phixk  = %.10f\n', phixk);
fprintf('  phixkp = %.10f\n', phixkp);

% -----------------------------------------------------------------------
% 7. Assemble A (6x6) and C (4x6)  (mleqadj.m lines 221-232)
% -----------------------------------------------------------------------
A_ss_mat = [ gammak,              gamma',     gamma0;
             zeros(4,1), P_mat,              P0_vec;
             0,      0,      0,      0,      0,      1];

B_mat = [ 0,0,0,0;
          Q_mat;
          0,0,0,0];

C_mat = [ [phiyk, phiyz, phiyl,    0, phiyg]+phiykp*Gamma(1:5)';
          [phixk,     0,     0,    0,     0]+phixkp*Gamma(1:5)';
          [philk, philz, phill,    0, philg]+philkp*Gamma(1:5)';
           0,0,0,0,1];

phi0  = Y0-C_mat*X0(1:5);
C_mat = [C_mat, phi0];

fprintf('\nA matrix (6x6):\n');
disp(A_ss_mat);
fprintf('C matrix (4x6)  [rows: y x l g | cols: log(k) log(z) taul taux log(g) const]\n');
disp(C_mat);
fprintf('phi0 = [%.10f, %.10f, %.10f, %.10f]\n', phi0(1), phi0(2), phi0(3), phi0(4));

% Quick sanity checks
fprintf('\nSanity checks:\n');
fprintf('  A(6,:) = ['); fprintf(' %.4f', A_ss_mat(6,:)); fprintf(' ]  (expected [0 0 0 0 0 1])\n');
fprintf('  A(2:5,2:5) == P_mat:  max_diff = %.2e\n', max(max(abs(A_ss_mat(2:5,2:5)-P_mat))));
fprintf('  A(1,1)   = gammak    = %.8f\n', A_ss_mat(1,1));
fprintf('  C(4,:)   = ['); fprintf(' %.4f', C_mat(4,:)); fprintf(' ]  (expected [0 0 0 0 1 0])\n');
fprintf('  Sbar(2)  = %.6f  (Table 8 target ~0.1336)\n', Sbar(1));
fprintf('  P(2,2)   = %.6f  (Table 8 target ~1.0011)\n', P_mat(2,2));

% -----------------------------------------------------------------------
% 8. Write CSVs
% -----------------------------------------------------------------------
dlmwrite([outdir '/A_bckm.csv'],        A_ss_mat,    'precision', '%.16e');
dlmwrite([outdir '/C_bckm.csv'],        C_mat,       'precision', '%.16e');
dlmwrite([outdir '/phi0_bckm.csv'],     phi0,        'precision', '%.16e');
dlmwrite([outdir '/gamma0_bckm.csv'],   gamma0,      'precision', '%.16e');
dlmwrite([outdir '/gamma_bckm.csv'],    [gammak; gamma]', 'precision', '%.16e');
dlmwrite([outdir '/Sbar_bckm.csv'],     Sbar,        'precision', '%.16e');
dlmwrite([outdir '/P_bckm.csv'],        P_mat,       'precision', '%.16e');
dlmwrite([outdir '/Qchol_bckm.csv'],    Q_mat,       'precision', '%.16e');

ss_row = [ks, ys, ls, xs, gs, cs, zs, tauls, tauxs, kls];
dlmwrite([outdir '/ss_calibrated.csv'], ss_row,      'precision', '%.16e');

% Partials table (one per row: name (as index), value)
partials = [philh; philk; philz; phill; philg; philkp;
            phiyk; phiyz; phiyl; phiyg; phiykp;
            phixk; phixkp;
            gammak; gamma0; gamma(1); gamma(2); gamma(3); gamma(4)];
dlmwrite([outdir '/partials_bckm.csv'], partials,    'precision', '%.16e');

fid = fopen([outdir '/partials_bckm_names.txt'], 'w');
names = {'philh','philk','philz','phill','philg','philkp', ...
         'phiyk','phiyz','phiyl','phiyg','phiykp', ...
         'phixk','phixkp', ...
         'gammak','gamma0','gamma_z','gamma_taul','gamma_taux','gamma_g'};
for i = 1:length(names)
  fprintf(fid, '%s,%.16e\n', names{i}, partials(i));
end
fclose(fid);

fprintf('\nCSVs written to %s/:\n', outdir);
fprintf('  A_bckm.csv              (6x6)\n');
fprintf('  C_bckm.csv              (4x6)\n');
fprintf('  phi0_bckm.csv           (4x1)\n');
fprintf('  gamma0_bckm.csv         (scalar)\n');
fprintf('  gamma_bckm.csv          (1x5: gammak gz gl gx gg)\n');
fprintf('  Sbar_bckm.csv           (4x1)\n');
fprintf('  P_bckm.csv              (4x4)\n');
fprintf('  Qchol_bckm.csv          (4x4 lower-tri Cholesky)\n');
fprintf('  ss_calibrated.csv       (ks ys ls xs gs cs zs tauls tauxs kls)\n');
fprintf('  partials_bckm.csv       (19x1)\n');
fprintf('  partials_bckm_names.txt (name,value)\n');
fprintf('\ndone.\n');
