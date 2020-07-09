%% Script to test different combinations of DoA and pseudorange residuals 
%% for spoof detection

clear variables;
close all;
rng(1);

saveFigures = false;

% three options: "or", "and", "joint"

% set max false alert probability
P_FAmax = 1e-7;

% set number of satellites
N = 12;
Nmin = 4; % minimum number of satellites considered

% set number of simulations
K = 1e5;

% set geometry (fixed)
az = linspace(0, 2*pi, N)';
el = linspace(pi/15, pi/2, N)';

% set pseudorange standard deviation, weight matrix
sigPr = 3; % m
w_pr = (1-0.5*cos(el)) ./ sigPr; % weight vector
% W = 1 ./ w.^2;

% set desired state fault / bias
xoffset = [0; 0; 10; 0];

% set DoA standard deviation
sigDoA = pi/6 * ones(N, 1);

% generate DoA measurements
yDoA_H0 = az + randn(N, K) .* sigDoA;
yDoA_H1 = randn(N, K) .* sigDoA;

% generate nominal pseudorange update measurements
y_H0 = diag(1./w_pr) * randn(N, K);

%% prepare doa, prr monitors

prr = LRT.pseudorangeResiduals(az, el, w_pr);
doa = LRT.SS_periodic(2*pi, az, zeros(N, 1), sigDoA, ones(N, 1));

%% update P_FAmax for number of hypothesis
% precalculate total number of considered subsets for threshold setting
numSubsets = 0;
for Nspo = Nmin:N
    subsets = binaryPermutations(N, Nspo)';
    numSubsets = numSubsets + size(subsets, 2);
end

doa_P_FAmax = [1/numSubsets, 1/sqrt(P_FAmax)/numSubsets, 1/(numSubsets+1)]*P_FAmax;
prr_P_FAmax = [1, 1/sqrt(P_FAmax), 1/(numSubsets+1)]*P_FAmax;
%% set prr thresholds

prr_thresholds = num2cell(prr.threshold(prr_P_FAmax));
[gamma_prr, gammaII_and, gammaII_or] = deal(prr_thresholds{:});

%% prepare and launch simulation


% loop over number of spoofed satellites
[P_MDmean, P_MDstd, P_MDmin, P_MDmax] = deal(zeros(N-Nmin+1, 5));
[eP_MDmean, eP_MDstd, eP_MDmin, eP_MDmax] = deal(zeros(N-Nmin+1, 5));
nSubsets = zeros(N-Nmin+1, 1);

for Nspo = Nmin:N

%% run empirical analysis of subsets
subsets = binaryPermutations(N, Nspo)';
nSubsets(Nspo-Nmin+1) = size(subsets, 2);
% initialize variables
[logLambdaI_H0, logLambdaI_H1, logLambdaII_H0, logLambdaII_H1] = ...
    deal(zeros(nSubsets(Nspo-Nmin+1), K));
[gammaI_or, gammaI_and, gamma_joint, gamma_doa, ...
    power_or, power_and, power_joint, power_doa, power_prr] = ...
    deal(zeros(nSubsets(Nspo-Nmin+1), 1));

fprintf('\n%i of %i satellites spoofed, %i subsets.\n', ...
    Nspo, N, nSubsets(Nspo-Nmin+1))

%% run for every possible subset
for ssI = 1:nSubsets(Nspo-Nmin+1)

Ispo = subsets(:, ssI);

if rem(ssI, 10)==0
    fprintf('%i ', ssI)
end
%% preprocessing calculations (pseudorange and DoA)

% calculate resulting pseudorange biases
y_b = prr.prBias(xoffset, Ispo);
% calculate noncentrality parameter
delta = prr.chi2statistic(y_b);

% create spoofed pseudorange updates
y_H1 = y_H0 + y_b;

% update doa object
doa.consideredSats = Ispo;

%% calculate thresholds and power

% calculate doa thresholds
doa_gammas = num2cell(doa.threshold(doa_P_FAmax));
[gamma_doa(ssI), gammaI_and(ssI), gammaI_or(ssI)] = deal(doa_gammas{:});

% create joint GLRT object
joint = LRT.GLRTcombined(doa, prr);

% calculate joint threshold
gamma_joint(ssI) = joint.threshold(P_FAmax / numSubsets);

% calcualte theoretical power
power_doa(ssI)   = doa.power(doa_P_FAmax(1));
power_prr(ssI)   = prr.power(prr_P_FAmax(1), delta);
power_and(ssI)   = doa.power(doa_P_FAmax(2)) ...
                 * prr.power(prr_P_FAmax(2), delta);
power_or(ssI)    = 1 - (1-doa.power(doa_P_FAmax(3))) ...
                    * (1-prr.power(prr_P_FAmax(3), delta));
power_joint(ssI) = joint.power(P_FAmax / numSubsets, delta);


%% calculate empiric metrics

% calculate log Lambda I | H_0
logLambdaI_H0(ssI, :) = doa.getLogLambda(yDoA_H0(Ispo, :));

% calculate log Lambda I | H_1
logLambdaI_H1(ssI, :) = doa.getLogLambda(yDoA_H1(Ispo, :));

% calculate log Lambda II | H_0
logLambdaII_H0(ssI, :) = prr.getLogLambda(y_H0);

% calculate log Lambda II | H_1
logLambdaII_H1(ssI, :) = prr.getLogLambda(y_H1);


end
fprintf('\n')
%% Calculate detection statistics

% empirically
detected_or = sum(logLambdaI_H1 < gammaI_or | logLambdaII_H1 < gammaII_or, 2);
detected_and = sum(logLambdaI_H1 < gammaI_and & logLambdaII_H1 < gammaII_and, 2);
detected_joint = sum(logLambdaI_H1 + logLambdaII_H1 < gamma_joint, 2);
detected_doa = sum(logLambdaI_H1 < gamma_doa, 2);
detected_prr = sum(logLambdaII_H1 < gamma_prr, 2);

falseAlert_or = sum(logLambdaI_H0 < gammaI_or | logLambdaII_H0 < gammaII_or, 2);
falseAlert_and = sum(logLambdaI_H0 < gammaI_and & logLambdaII_H0 < gammaII_and, 2);
falseAlert_joint = sum(logLambdaI_H0 + logLambdaII_H0 < gamma_joint, 2);
falseAlert_doa = sum(logLambdaI_H0 < gamma_doa, 2);
falseAlert_prr = sum(logLambdaII_H0 < gamma_prr, 2);

eP_MD = 1-[detected_doa, detected_prr, detected_and, detected_or, detected_joint]/K;
eP_MDmean(Nspo-Nmin+1, :) = mean(eP_MD, 1);
eP_MDstd(Nspo-Nmin+1, :) = std(eP_MD, 0, 1);
eP_MDmin(Nspo-Nmin+1, :) = min(eP_MD, [], 1);
eP_MDmax(Nspo-Nmin+1, :) = max(eP_MD, [], 1);

% analytically
P_MD = 1-[power_doa, power_prr, power_and, power_or, power_joint];
P_MDmean(Nspo-Nmin+1, :) = mean(P_MD, 1);
P_MDstd(Nspo-Nmin+1, :) = std(P_MD, 0, 1);
P_MDmin(Nspo-Nmin+1, :) = min(P_MD, [], 1);
P_MDmax(Nspo-Nmin+1, :) = max(P_MD, [], 1);

eP_FA = [falseAlert_doa, falseAlert_prr, ...
    falseAlert_and, falseAlert_or, falseAlert_joint]/K;

fprintf(['empirical P_MD: %.2f %% (doa); %.2f %% (prr) %.2f %% (and) ', ...
    '%.2f %% (or) %.2f %% (joint)\n'], ...
    100*eP_MDmean(Nspo-Nmin+1, :))
fprintf(['theoretical P_MD: %.2f %% (doa); %.2f %% (prr) %.2f %% (and) ', ...
    '%.2f %% (or) %.2f %% (joint)\n'], ...
    100*P_MDmean(Nspo-Nmin+1, :))
fprintf(['empirical P_FA (*P_FAmax): %.2f (doa); %.2f (prr) %.2f (and) ', ...
    '%.2f (or) %.2f (joint)\n'], ...
    mean(eP_FA, 1) / P_FAmax)

if Nspo == 7
    %% plot single case
    [~, ssI] = max(eP_MD(:, 2)); % most difficult for prr
    fs = 18;
    f1 = figure; hold on; grid on;
    plot(logLambdaI_H0(ssI, :), logLambdaII_H0(ssI, :), '+')
    plot(logLambdaI_H1(ssI, :), logLambdaII_H1(ssI, :), 'o')
    ylim([max(min([logLambdaII_H1(ssI, :), -50]), -100), 0])
    xlim([min(logLambdaI_H1(ssI, :)), max(logLambdaI_H0(ssI, :))])
    % plot "or" threshold
    plot([gammaI_or(ssI)*ones(1, 2), f1.CurrentAxes.XLim(2)], ...
        [f1.CurrentAxes.YLim(2), gammaII_or*ones(1, 2)], ...
        'k', 'LineWidth', 1.3)
    % plot "and" threshold
    plot([f1.CurrentAxes.XLim(1), gammaI_and(ssI)*ones(1, 2)], ...
        [gammaII_and*ones(1, 2), f1.CurrentAxes.YLim(1)], ...
        'k--', 'LineWidth', 1.3)
    % plot "joint" threshold
    plot([gamma_joint(ssI), f1.CurrentAxes.XLim(2)], ...
        [0, gamma_joint(ssI) - f1.CurrentAxes.XLim(2)], ...
        'k-.', 'LineWidth', 1.5)
    legend('$H_0$', '$H_1$', '$\gamma_{OR}$', '$\gamma_{AND}$', ...
        '$\gamma_{joint}$', 'Location', 'best', ...
        'FontSize', fs, 'Interpreter', 'latex')
    xlabel('$\log\Lambda(y_{DoA})$', ...
        'FontSize', fs, 'Interpreter', 'latex')
    ylabel('$\log\Lambda(y_{prr})$', ...
        'FontSize', fs, 'Interpreter', 'latex')
    title([num2str(Nspo), ' out of ', num2str(N), ' satellites spoofed'], ...
        'FontSize', fs, 'Interpreter', 'latex')
end
end
%% Plot results
fs = 18;

% show course of P_MD over detections
fAll = figure; hold on; grid on;
plot(Nmin:N, 100*P_MDmean(:, 1:2)', '--', 'LineWidth', 1.5)
plot(Nmin:N, 100*P_MDmean(:, 3:5)', 'LineWidth', 1.5)
legend('DoA LRT', 'Residuals $\chi^2$ test', '"and"', '"or"', '"joint"', ...
    'Location', 'best', 'FontSize', fs, 'Interpreter', 'latex')
ylim([0, 100])
fAll.CurrentAxes.FontSize = fs-4;
xlabel(['Number of spoofed satellites (out of ', num2str(N), ')'], ...
    'FontSize', fs, 'Interpreter', 'latex')
ylabel('$E_{s\in S\sim U} [P_{MD}]$ in $\%$', ...
    'FontSize', fs, 'Interpreter', 'latex')

fMax = figure; hold on; grid on;
plot(Nmin:N, 100*P_MDmax', ...
    'LineWidth', 1.5)
ylim([0, 100])
fMax.CurrentAxes.FontSize = fs-4;
legend('DoA LRT', 'Residuals $\chi^2$ test', '"and"', '"or"', '"joint"', ...
    'Location', 'best', 'FontSize', fs, 'Interpreter', 'latex')
xlabel(['Number of spoofed satellites (out of ', num2str(N), ')'], ...
    'FontSize', fs, 'Interpreter', 'latex')
ylabel('$\max_{s \in S}$ $P_{MD}$ in $\%$', ...
    'FontSize', fs, 'Interpreter', 'latex')

fMaxEmp = figure; hold on; grid on;
plot(Nmin:N, 100*eP_MDmax', ...
    'LineWidth', 1.5)
ylim([0, 100])
fMaxEmp.CurrentAxes.FontSize = fs-4;
legend('DoA LRT', 'Residuals $\chi^2$ test', '"and"', '"or"', '"joint"', ...
    'Location', 'best', 'FontSize', fs, 'Interpreter', 'latex')
xlabel(['Number of spoofed satellites (out of ', num2str(N), ')'], ...
    'FontSize', fs, 'Interpreter', 'latex')
ylabel('$\max_{s \in S}$ $eP_{MD}$ in $\%$', ...
    'FontSize', fs, 'Interpreter', 'latex')

%% 
fSubplot = figure;
subplot(2, 1, 1); hold on; grid on;
plot(Nmin:N, 100*eP_MDmax', ...
    'LineWidth', 1.5)
ylim([0, 100])
fSubplot.CurrentAxes.FontSize = fs-4;
ylabel('$\max_{s \in S}$ $eP_{MD}$ in $\%$', ...
    'FontSize', fs, 'Interpreter', 'latex')

subplot(2, 1, 2); hold on; grid on;
plot(Nmin:N, eP_MDmax', ...
    'LineWidth', 1.5)
fSubplot.CurrentAxes.FontSize = fs-4;
fSubplot.CurrentAxes.YScale = 'log';
legend('DoA LRT', 'Residuals $\chi^2$ test', 'AND', 'OR', '"joint"', ...
    'Location', 'south', 'FontSize', fs, 'Interpreter', 'latex')
xlabel(['Number of spoofed satellites (out of ', num2str(N), ')'], ...
    'FontSize', fs, 'Interpreter', 'latex')
ylabel('$\max_{s \in S}$ $eP_{MD}$', ...
    'FontSize', fs, 'Interpreter', 'latex')
fSubplot.Position(4) = 600;

%% Print maximum integrity risk for each detection method:
fprintf(['Maximum integrity risk (P_MD): ', ...
    '%.2f %% (doa); %.2f %% (prr) %.2f %% (and) ', ...
    '%.2f %% (or) %.2f %% (joint)\n'], ...
    100*max(P_MDmax, [], 1))

disp(['Reduction of max P_MD of joint vs. OR: ', ...
    num2str(round((1 - max(eP_MDmax(:, end)) ./ max(eP_MDmax(:, 4)))*100, 2)), ...
    '%'])
