%%Simulation study combining power and distortion measurements for spoofing
%%detection
%   
%   Computes the detection power for detection tests using the power and
%   distortion metrics. Follows algorithms described in
%   
%   Rothmaier, F., Chen, Y., Lo, S., & Walter, T. (2021). A Framework for 
%       GNSS Spoofing Detection through Combinations of Metrics. IEEE 
%       Transactions on Aerospace and Electronic Systems. 
%       https://doi.org/10.1109/TAES.2021.3082673
%   
%   The simulation is described and analyzed in
%   Rothmaier, F., Lo, S., Phelts, E., & Walter, T. (2021). GNSS Spoofing 
%       Detection through Metric Combinations: Calibration and Application 
%       of a general Framework. Proceedings of the 34th International 
%       Technical Meeting of the Satellite Division of the Institute of 
%       Navigation, ION GNSS+ 2021.
%   
%   Please cite the reference if using this material.
%   
%   This code and the repository it is in is provided free of charge under
%   an MIT license. For questions please contact Fabian Rothmaier
%   fabianr@stanford.edu.
%   
%   Written by Fabian Rothmaier, Stanford University, 2021


%% prepare workspace
clear variables;
close all;

rng(1); % for reproducability

%% set simulation paramters

CN0nom = 40; % nominal signal strength in dB-Hz
sigP = 2; % standard deviation of power metric in dB

ELdist = 0.052; % early/late correlator positions in chips from prompt

P_FA = 1e-8; % false alert probability

deltaCN0s = 1:20; % different spoofer power advantages in dB

N = 1e6; % number of monte carlo runs at every epoch

offsets = 0:0.01:(1+ELdist+0.01); % offsets of spoofer's triangle in code chips

dx = 0.001;
x = -1-min(offsets):dx:1+max(offsets); % code chip range around authentic peak
xTI = round(ELdist/dx) * [-1 0 1]; % correlation tab spacing in units of x

% calculate height of nominal signal peak at each epoch
nomPeak = expectedPeaksize(CN0nom);
nomCorrFun = nomPeak .* Rcorr(x'); % nominal correlation triangle

%% run simulation during lift-off
deltaMetric = zeros(length(deltaCN0s), length(offsets));
peakSize = zeros(length(deltaCN0s), length(offsets));

[meanLLD, maxLLD, meanLLP, maxLLP] = deal(zeros(1, length(deltaCN0s)));

for deltaI = 1:length(deltaCN0s)
    deltaCN0 = deltaCN0s(deltaI); % spoofer power advantage


    %% set peak sizes
    CN0spo = CN0nom + deltaCN0;

    spoPeak = expectedPeaksize(CN0spo);
    
    % compute spoofed correlation function
    spoCorrFun = spoPeak .* Rcorr(x' - offsets);
    
    % sum with nominal values
    sumCorrFun = nomCorrFun  + spoCorrFun;

    %% compute nominal peak, spoofed peak
    
    peakIs = zeros(size(offsets));
    
    deltaMetric = zeros(N, length(offsets));
    measCN0 = zeros(N, length(offsets));
    
    for ii = 1:length(offsets)

        % find peak, assume perfect tracking
        [peakSize(deltaI, ii), peakIs(ii)] = max(sumCorrFun(:, ii));

        % run a Monte Carlo for Delta metric measurement
        corrTapStd = deltaSigmaVal(expectedCN0(peakSize(deltaI, ii)))/sqrt(2);
        corrTapNoise = corrTapStd*randn(N, length(xTI)); % is normalized by peaksize
        corrTapVals = sumCorrFun(peakIs(ii)+xTI, ii)' + peakSize(deltaI, ii)* corrTapNoise;
        
        deltaMetric(:, ii) = (corrTapVals(:, 1) - corrTapVals(:, 3)) ./ corrTapVals(:, 2);
        
        % get measured CN0 (mostly unaffected by noise on correlator taps)
        measCN0(:, ii) = expectedCN0(corrTapVals(:, 2));
        
    end
    
    % calculate log Lambda of the delta metric
    sigDelta = deltaSigmaVal(measCN0);
    logLaD = (deltaMetric ./ sigDelta).^2;
    
    % calculate log Lambda of the power metric, now do Monte Carlo for
    % noise on power metric (excluded before for comp. efficiency, does not
    % affect delta metric)
    % This includes extra noise on power measurement.
    logLaP = ((measCN0 + sigP*randn(N, length(offsets)) - CN0nom) / sigP).^2;
    
    % calculate alarms
    alarm.D = sum(logLaD > chi2inv(1-P_FA, 1), 1);
    alarm.P = sum(logLaP > chi2inv(1-P_FA, 1), 1);
    alarm.Dor = sum(logLaD > chi2inv(1-P_FA/2, 1), 1);
    alarm.Por = sum(logLaP > chi2inv(1-P_FA/2, 1), 1);
    alarm.OR = sum(logLaD > chi2inv(1-P_FA/2, 1) | logLaP > chi2inv(1-P_FA/2, 1), 1);
    alarm.Dand = sum(logLaD > chi2inv(1-sqrt(P_FA), 1), 1);
    alarm.Pand = sum(logLaP > chi2inv(1-sqrt(P_FA), 1), 1);
    alarm.AND = sum(logLaD > chi2inv(1-sqrt(P_FA), 1) & logLaP > chi2inv(1-sqrt(P_FA), 1), 1);
    alarm.joint = sum(logLaP + logLaD > chi2inv(1-P_FA, 2), 1);
    
    for metric = fields(alarm)'
        % minimum P_MD at best epoch
        P_MDmin.(metric{1})(deltaI) = 1 - max(alarm.(metric{1})) / N;
        % average P_MD across all epochs
        P_MDmean.(metric{1})(deltaI) = 1 - mean(alarm.(metric{1})) / N;
    end

end

%% plot results
plt = PlotLatexStyle;

% compute detection with OR, AND, joint metric

% plot number of alarms
figure; hold on; grid on;
plot(deltaCN0s, N*(1-P_MDmean.P), '-+', 'LineWidth', 1.2)
plot(deltaCN0s, N*(1-P_MDmean.D), 'LineWidth', 1.2)
plot(deltaCN0s, N*(1-P_MDmean.AND), '--', 'LineWidth', 1.2)
plot(deltaCN0s, N*(1-P_MDmean.OR), '-.', 'LineWidth', 1.2)
plot(deltaCN0s, N*(1-P_MDmean.joint), '-d', 'LineWidth', 1.2)

plt.latexLegend('P', 'D', 'AND', 'OR', 'joint')

xlabel('Spoofer power advantage in (dB)', plt.axisLabelArgs{:})
ylabel('Number of alarms', plt.axisLabelArgs{:})


% plot average P_MD in %
figure; hold on; grid on;
plot(deltaCN0s, P_MDmean.P*100, '-+', 'LineWidth', 1.2)
plot(deltaCN0s, P_MDmean.D*100, 'LineWidth', 1.2)
plot(deltaCN0s, P_MDmean.AND*100, '--', 'LineWidth', 1.2)
plot(deltaCN0s, P_MDmean.OR*100, '-.', 'LineWidth', 1.2)
plot(deltaCN0s, P_MDmean.joint*100, '-d', 'LineWidth', 1.2)

plt.latexLegend('P', 'D', 'AND', 'OR', 'joint')

xlabel('Spoofer power advantage in (dB)', plt.axisLabelArgs{:})
ylabel('Average $P_{MD}$ during lift-off in $\%$ ', plt.axisLabelArgs{:})


% Plot P_MD on logscale
fLog = figure; hold on; grid on;

plot(deltaCN0s, P_MDmean.P, '-+', 'LineWidth', 1.2)
plot(deltaCN0s, P_MDmean.D, 'LineWidth', 1.2)
plot(deltaCN0s, P_MDmean.AND, '--', 'LineWidth', 1.2)
plot(deltaCN0s, P_MDmean.OR, '-.', 'LineWidth', 1.2)
plot(deltaCN0s, P_MDmean.joint, '-d', 'LineWidth', 1.2)
fLog.CurrentAxes.FontSize = plt.fs-4;
fLog.CurrentAxes.YScale = 'log';

plt.latexLegend('P', 'D', 'AND', 'OR', 'joint')
xlabel('Spoofer power advantage in (dB)', plt.axisLabelArgs{:})
ylabel('Average $P_{MD}$ during lift-off', plt.axisLabelArgs{:})


%% helper functions
function ps = expectedPeaksize(CN0, N0)
%Calculate the theoretical peak size based on CN0, N0.

if nargin < 2
    N0 = -203;
end

ps = sqrt(10.^((CN0 + N0)/10))*10^12 .* 500;

end


function eCN0 = expectedCN0(peaksize, N0)
%Compute the expected CN0 as a function of correlation peak size.

if nargin < 2
    N0 = -203;
end

% compute the expected C/N0
eCN0 = 10 * log10((2*peaksize /1e15).^2) - N0;

end