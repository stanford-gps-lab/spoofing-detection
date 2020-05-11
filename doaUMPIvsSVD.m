%% DoA UMPI hypothesis test vs. svd based detection comparison
% 
% This script is to compare the performance of the NP-test based approach
% to the SVD based approach to detect spoofing.
% This is an application example of the UMPI hypothesis test presented in 
% 
% Rothmaier, F. et al. GNSS Spoof Detection through Spatial Processing.
% NAVIGATION, 2020
% 
% The SVD based approach is presented in:
% Appel, M. et al. Experimental validation of GNSS repeater detection 
% based on antenna arrays for maritime applications. CEAS Sp. J. 11, 
% 7–19 (2019).
% 
% In this script we generate random azimuth and elevation measurements and
% perform the NP-test based detection. We plot the detection metric over
% time to display the detection power.
% 
% 
% Written by Fabian Rothmaier, 2020
% 
% Published under MIT license.

%% prepare work environment

% clear and prepare workspace
close all;
clear variables;

% seed random number generator for repeatable results
rng(2);
fs = 20; % figure font size

%% User inputs

% number of measurement epochs
K = 2e3;

% number of satellites in view
N = 10;

% choose scenario (1 => 15 deg cov, 2 => 25 deg cov), see below
scenario = 1;

% randomized measurement uncertainty for each DoA?
randomSigma = true; % if false it's governed by covDeg variable

% set false alert probability 10^(x)
pFAmaxPot = -8;

%% define scenario
if scenario == 1
    linestyle = 'g';
    spoofedEpochs = [200:400, 650:800, 1000:1500];
    covDeg = 15;
else
    linestyle = 'b+-';
    spoofedEpochs = [200:400, 650:800, 1000:2000];
    covDeg = 25;
end

% update K to include spoofed Epochs
K = max(max(spoofedEpochs), K);

% create logical indices
spoEps = false(1, K);
spoEps(spoofedEpochs) = true;

% number of sequential measurements to consider
nSeq = 3;

%% describe constellation for simulation

azSpoof = 2; % azimuth direction of spoofer
elSpoof = 0.8; % elevation of spoofer

% azTrue = repmat(2 * pi * rand(N, 1), 1, K); % uniform random between 0 and 2 pi
% elRand = repmat(pi / 2 * rand(N, 1), 1, K); % uniform random between 0 and pi
% azTrue = repmat(2 * pi * [0; 0.02], 1, K); % specific between 0 and 2 pi
% elRand = repmat(pi / 2 * [0.5; 0.5], 1, K); % specific between 0 and pi
azTrue = 2*pi * rand(N, K); % random at each epoch
elRand = pi/2 * rand(N, K); % random at each epoch
elTrue = max(min(elRand, (pi/2-1e-4)*ones(N, K)), zeros(N, K));

% define ephemeris-based DoA unit vectors
ephemerisUVs = azel2enu(azTrue, elTrue);

% generate measurement noise
if ~randomSigma
    sigma = repmat(sqrt(covDeg)*pi/180, N, K);
else
    sigma = repmat(sqrt(15:14+N)'*pi/180, 1, K); % different sigma_n
end

% add spoofer in desired time slots
azTrueS = azTrue; azTrueS(:, spoEps) = azSpoof;
elTrueS = elTrue; elTrueS(:, spoEps) = elSpoof;

%% add noise to generate measurements

[azMeas, elMeas] = doa2Dnoise(azTrueS, elTrueS, sigma);

% convert to DoAs unit vectors
measuredUVs = azel2enu(azMeas, elMeas);

%% prepare simulation
% precompute quartile value
PhiT = norminv(10^pFAmaxPot);

% preallocate result variables
[SigCA, logLambdaCAoffset, logLambdaNormCAseq] = deal(zeros(1, K));
Q = zeros(3, K);
DOP = zeros(4, K);

% normalized arc errors
normErrors = zeros((2*N - 3), K);
conditioning = zeros(2*N - 3, 1);

S_barCA = zeros(2*N - 3, 2*N - 3, K);

%% run simulation
for k = 1:K
    
    % rotate antenna arbitrarily (does not change the result)
%     att = [(rand(1, 2)-0.5)*pi, rand(1)*2*pi];
%     R1 = [1 0 0; ...
%           0 cos(att(1)) sin(att(1)); ...
%           0 -sin(att(1)) cos(att(1))];
%     R2 = [cos(att(2)) 0 -sin(att(2)); ...
%           0 1 0; ...
%           sin(att(2)) 0 cos(att(2))];
%     R3 = [cos(att(3)) sin(att(3)) 0; ...
%           -sin(att(3)) cos(att(3)) 0; ...
%           0 0 1];
    
    % calculate normalized singular values
    Q(:, k) = svd(measuredUVs(:, :, k)*ephemerisUVs(:, :, k)') / N;
     
    % select arcs and define measurement covariance
    GCAeph = GreatCircleArcSelection(ephemerisUVs(:, :, k), sigma(:, k));
    [arcSel, S_barCA(:, :, k), conditioning(k)] = GCAeph.findArcs(100);
    phi_bar = GCAeph.getGCAs(arcSel);
    [U, S, V] = svd(S_barCA(:, :, k));
    % get measured values
    y_bar = GCAeph.computeGCAs(arcSel, measuredUVs(:, :, k))';
    
    SigCA(k) = phi_bar' * (V/S*U' * phi_bar);
    normErrors(:, k) = (sqrt(S) \ U') * (y_bar - phi_bar);
    % distance of logLambda from expected mean
    logLambdaCAoffset(k) = phi_bar' * (V/S*U' * (y_bar - phi_bar));
    
    % sequential Lambda (consider nSeq epochs)
    seqEp = max([1, k-nSeq+1]) : k;
    logLambdaNormCAseq(k) = sum(logLambdaCAoffset(seqEp)) ...
        / sqrt(sum(SigCA(seqEp)));
    
    % store geometry in form of DOP
    G = [-sin(azTrue(:, k)).*cos(elTrue(:, k)), ...
         -cos(azTrue(:, k)).*cos(elTrue(:, k)), ...
         -sin(elTrue(:, k)), ones(N, 1)];
    DOP(:, k) = diag(inv(G'*G));

end

disp(['Average P_MD = ', num2str(mean(normcdf(sqrt(SigCA) + PhiT, 'upper')))])

%% Plot results

% histogram of nominal log Lambdas
figure; hold on; grid on;
histogram(normErrors(:, ~spoEps), round(K/50), ...
    'normalization', 'pdf', ...
    'FaceColor', 'g', 'EdgeColor', 'none');
plot(linspace(PhiT, -PhiT), normpdf(linspace(PhiT, -PhiT)), ...
    'Color', [0 0.5 0], 'LineWidth', 1.5)
xlabel('$\bar{y} - \bar{\phi}$ normalized', ...
    'FontSize', fs, 'Interpreter', 'latex')
ylabel('pdf', 'FontSize', fs, 'Interpreter', 'latex')

% histogram of nominal log Lambdas using central angle approach
fH = figure; hold on; grid on;
histogram(logLambdaCAoffset(~spoEps)./sqrt(SigCA(~spoEps)), ...
    min(round(K/50), 200), 'normalization', 'pdf', ...
    'FaceColor', 'g', 'EdgeColor', 'none');
plot(linspace(PhiT, -PhiT), normpdf(linspace(PhiT, -PhiT)), ...
    'Color', [0 0.5 0], 'LineWidth', 1.5)
line(PhiT*ones(1, 2), fH.CurrentAxes.YLim, 'Color', 'k', 'LineWidth', 2)
xlabel('$\log\Lambda(\bar{y}) | H_0$ normalized', ...
    'FontSize', fs, 'Interpreter', 'latex')
ylabel('pdf', 'FontSize', fs, 'Interpreter', 'latex')

% skyplot of ephemeris based and measured data
figure;
plot_skyplot(azTrue(:, 1), elTrue(:, 1), azMeas, elMeas)

% snapshot and sequential normalized log Lambdas using great circle arcs
fig1 = figure; hold on; grid on;
plot(logLambdaCAoffset ./ sqrt(SigCA), linestyle, 'LineWidth', 1.5);
plot(fig1.CurrentAxes.YLim(2)*(logLambdaCAoffset ./ sqrt(SigCA) < PhiT) ...
    + fig1.CurrentAxes.YLim(1)*(logLambdaCAoffset ./ sqrt(SigCA) >= PhiT), ...
        'b', 'LineWidth', 1.5)
fig1.CurrentAxes.Children = flipud(fig1.CurrentAxes.Children);
line([0 K], [PhiT PhiT], 'Color', 'k') % decision threshold
text(K*0.75, PhiT-2, ['$P_{FA_{max}} = 10^{', num2str(pFAmaxPot), '}$'], ...
    'FontSize', fs-6, 'Interpreter', 'latex')
title(['$\sigma^2 = ', num2str(covDeg), '$ deg$^2$, ', ...
    'snapshot-based'], ...
    'FontSize', fs, 'Interpreter', 'latex')
xlabel('Simulation run', 'FontSize', fs, 'Interpreter', 'latex')
ylabel('Normalized $\log\Lambda(\bar{y}_t)$', ...
    'FontSize', fs, 'Interpreter', 'latex')

% plot sequential log Lambda for comparison
fig2 = figure; hold on; grid on;
p1 = plot(logLambdaCAoffset ./ sqrt(SigCA), [linestyle(1), ':']);
p2 = plot(logLambdaNormCAseq, linestyle, 'LineWidth', 1.5);
ylim([-50 10])
plot(fig2.CurrentAxes.YLim(2)*(logLambdaNormCAseq < PhiT) ...
    + fig2.CurrentAxes.YLim(1)*(logLambdaNormCAseq >= PhiT), ...
        'b', 'LineWidth', 1.5)
fig2.CurrentAxes.Children = flipud(fig2.CurrentAxes.Children);
line([0 K], [PhiT PhiT], 'Color', 'k') % decision threshold
text(K*0.75, PhiT-2, ['$P_{FA_{max}} = 10^{', num2str(pFAmaxPot), '}$'], ...
    'FontSize', fs-6, 'Interpreter', 'latex')
legend([p1; p2], {'snapshot'; 'sequential'}, ...
    'Location', 'southeast', 'FontSize', fs-8, 'Interpreter', 'latex')
title(['$\sigma^2 = ', num2str(covDeg), '$ deg$^2$, ', ...
    '$', num2str(nSeq), '$ sequential epochs'], ...
    'FontSize', fs, 'Interpreter', 'latex')
xlabel('Simulation run', 'FontSize', fs, 'Interpreter', 'latex')
ylabel(['Normalized $\log\Lambda(\bar{y}_{(t-', num2str(nSeq-1), '):t})$'], ...
    'FontSize', fs, 'Interpreter', 'latex')

% to compare: normalized sum of singular values
figure; hold on; grid on;
plot(sum(Q), linestyle, 'LineWidth', 1.5)
ylim([0 1.1])
xlabel('Simulation run', 'FontSize', fs, 'Interpreter', 'latex')
ylabel('Normalized $\sum\sigma_i$', ...
    'FontSize', fs, 'Interpreter', 'latex')


