%% Missed Detection probability of GLRT and UMPI test
% 
% Script to examine the P_MD of the GLRT and UMPI tests when run on the
% same Direction of Arrival (DoA) measurements. The script runs through all
% possble subsets of spoofed satellites between 4 and all satellites in the
% defined constellation and calculates the missed detection probability of
% the simple vs. simple UMPI test and Composite vs. Composite GLRT.
% 
% Written by Fabian Rothmaier, 2020
% 
% Published under MIT license.

clear variables;
close all;
rng(2); % for reproducability

%% define constellation, decision threshold
% set min, max number of satellites
Nmax = 9;
Nmin = 4;

K = 1e6; % number of epochs

% false alert probability
P_FAmax = 1e-7;

% define satellite constellation
azTrue = 2 * pi * rand(Nmax, 1); % uniform random between 0 and 2 pi
elRand = pi / 2 * rand(Nmax, 1); % uniform random between 0 and pi
% azTrue = 2 * pi * [0; 0.02]; % specific between 0 and 2 pi
% elRand = pi / 2 * [0.5; 0.5]; % specific between 0 and pi/2

% avoid singularity at 90deg
elTrue = min(elRand, (pi/2-1e-4)*ones(Nmax, 1));

% set DoA standard deviation
sigDoA = 12/180*pi * ones(Nmax, 1);

% define DoA unit vectors
ephemUVs = azel2enu(azTrue, elTrue);

%% create spoofed measurements

[azMeas, elMeas] = doa2Dnoise(zeros(Nmax, K), zeros(Nmax, K), sigDoA);

% convert to DoAs unit vectors
measuredUVs = azel2enu(azMeas, elMeas);

%% Prepare subset hypothesis
% count number of subsets totally considered
nHypothesis = 0;
for N = Nmin:Nmax
    nHypothesis = nHypothesis + size(binaryPermutations(Nmax, N), 1);
end

% calculate P_FA quantile per test
PhiT = norminv(P_FAmax / nHypothesis);
PhiT_chi2 = chi2inv(1 - P_FAmax/nHypothesis, 2*N-3);

%% run over all subsets

for N = Nmax:-1:Nmin
    
    subsets = binaryPermutations(Nmax, N)';
    nSubsets = size(subsets, 2);
    
    [SigCA, lambda, P_MD_CAempiric] = deal(zeros(nSubsets, 1));
    
    for ssI = 1:nSubsets
        % pick subset
        Ispo = subsets(:, ssI);

        %% compute central angles
        % (use GCA-Selection class)
        GCA = GreatCircleArcSelection(ephemUVs(:, Ispo), sigDoA(Ispo));
        [arcSel, S_barCA, minCond, ki] = GCA.findArcs(100);
        phi_bar = GCA.getGCAs(arcSel);
        SigCA(ssI) = phi_bar' * (S_barCA \ phi_bar);
                                      
        % get measured values
        y_bar = GCA.computeGCAs(arcSel, measuredUVs(:, Ispo, :));
        
        % run Monte Carlo to compare against theoretical value:
        % log Lambda - mu_LambdaH0:
        logLambdaCA = phi_bar' * (S_barCA \ (y_bar - phi_bar));
        % normalize by sigma_LambdaH0, count missed detections
        P_MD_CAempiric(ssI) = sum(logLambdaCA / sqrt(SigCA(ssI)) >= PhiT) / K;
        
        
        %% comparison: GLRT approach
        
%         % svd if all measurements are considered for attitude computation
%         B = unitVectors;
%         B(:, Ispo) = measuredUVs(:, Ispo, 1);
%         % calculate rotation matrix
%         [U,S,V] = svd(B * unitVectors');
        
        % calculate rotation matrix
        [U,S,V] = svd(measuredUVs(:, Ispo, 1) * ephemUVs(:, Ispo)');
        R = U * diag([1, 1, det(U*V')]) * V';
        
        lambda(ssI) = sum(( acos( ...
            diag(ephemUVs(:, Ispo)' * R' * measuredUVs(:, Ispo, 1))) ...
            ./ sigDoA(Ispo)).^2, 1);
    end
    
    % compute P_MD stats
    P_MD_CAvalues = normcdf(sqrt(SigCA) + PhiT, 'upper');
    P_MD_CA.max(N-Nmin+1) = max(P_MD_CAvalues);
    P_MD_CA.mean(N-Nmin+1) = mean(P_MD_CAvalues);
    P_MD_CA.min(N-Nmin+1) = min(P_MD_CAvalues);
    P_MD_CA.std(N-Nmin+1) = std(P_MD_CAvalues);
    P_MD_CA.emp(N-Nmin+1) = max(mean(P_MD_CAempiric), 1/K);
    
    P_MDchi2values = ncx2cdf(PhiT_chi2, 2*N-3, lambda);
    P_MDchi2.max(N-Nmin+1) = max(P_MDchi2values);
    P_MDchi2.mean(N-Nmin+1) = mean(P_MDchi2values);
    P_MDchi2.min(N-Nmin+1) = min(P_MDchi2values);
    P_MDchi2.std(N-Nmin+1) = std(P_MDchi2values);
    
end

%% plot results
fs = 18;

f1 = figure;
semilogy(Nmin:Nmax, P_MDchi2.mean, 'LineWidth', 2)
hold on; grid on;
semilogy(Nmin:Nmax, P_MD_CA.emp, '--', 'LineWidth', 2)

f1.CurrentAxes.FontSize = fs - 4;
xlabel('Number of spoofed satellites', ...
    'FontSize', fs, 'Interpreter', 'latex')
ylabel('Average $P_{MD}$', ...
    'FontSize', fs, 'Interpreter', 'latex')
legend('GLRT', 'UMPI (empirical)', ...
    'Location', 'best', 'FontSize', fs, 'Interpreter', 'latex')
f1.Position(4) = 250;
