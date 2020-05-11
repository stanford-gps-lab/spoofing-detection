classdef GreatCircleArcSelection
    %GreatCircleArcSelection Choose selection of arcs among DoAs
    %   Chooses a subset of arcs that results in a well conditioned
    %   covariance matrix and large power of the Likelihood Ratio Test 
    %   (LRT). Choosing the best selection of arcs is NP-hard. This class
    %   efficiently finds a reasonable good selection that results in a low
    %   conditioning number of the arc measurement covariance matrix and in
    %   a large Mahalanobis distance arcs' / covar * arcs. A larger
    %   Mahalanobis distance results in a more powerful LRT weather some
    %   measured arcs match the true arcs or are all zero.
    %   
    %   Properties are:
    %       doaUnitVectors  unit vectors of DoA measurements
    %       N               number of DoA measurements
    %       sigmas          DoA measurement standard deviations
    %       gcas            matrix of all possible great circle arcs
    %   
    %   Methods:
    %       obj = GreatCircleArcSelection(doas, sigmas)
    %           constructor
    %       [arcs, S_bar, minCond, ki] = findArcs(obj, K, maxCond)
    %           finds best combination of great circle arcs
    %       S_bar_sel = getCovarianceMatrix(obj, gcaSelection)
    %           construct Covariance matrix
    %       gcas = getGCAs(obj, gcaSelection)
    %           computes vector of great circle arcs
    %       cosTheta = sphericalCosine(obj, doaIndices)
    %           Computes spherical cosine between two great circle arcs
    % 
    %   Written by Fabian Rothmaier, 2020
    %   
    %   Published under MIT license
    
    properties
        doaUnitVectors  % unit vectors of DoA measurements
        N               % number of DoA measurements
        sigmas          % DoA measurement standard deviations
        gcas            % matrix of all great circle arcs
    end
    
    methods
        function obj = GreatCircleArcSelection(doas, sigmas)
            %GreatCircleArcSelection(doas, sigmas)
            %   Construct an instance to select great circle arcs.
            
            % make sure doas are column vectors
            [n, m] = size(doas);
            if n ~= 3 && m == 3
                doas = doas';
            elseif ~any([n, m] == 3)
                error('Passed doas need to be 3xN matrix.')
            end
            
            % store doa unit vectors
            obj.doaUnitVectors = doas ./ sqrt(sum(doas.^2, 1));
            % store number of doas
            obj.N = size(doas, 2);
            
            % store doa standard deviations
            obj.sigmas = sigmas;
            
            % store matrix of great circle arcs
            gcaOptions = nchoosek(1:obj.N, 2); % possible combinations
            % compute arcs
            gcaList = computeGCAs(obj, gcaOptions);
            % reshape to matrix
            obj.gcas = zeros(obj.N);
            obj.gcas(sub2ind(size(obj.gcas), gcaOptions(:, 1), gcaOptions(:, 2))) = ...
                gcaList;
            obj.gcas = obj.gcas + obj.gcas';
        end
        
        
        function [arcs, S_bar, minCond, ki] = findArcs(obj, K, maxCond)
            %obj.getCovariance(K) finds best combination of great circle
            %   arcs using K samples. A larger number of samples results in
            %   a better conditioned covariance matrix and larger
            %   Mahalanobis distance but causes longer computation time.
            %   
            %   Input:
            %   K           max number of samples (default 100)
            %   maxCond     max conditioning number of S_bar
            %   
            %   Output:
            %   arcs    2N-3 x 2 matrix of doa indices defining arcs
            %   S_bar   2N-3 x 2N-3 covariance matrix
            %   minCond conditioning of number of covariance matrix
            %   ki      number of samples drawn until success
            
            % random sample of possible arc combinations
            
            if nargin < 3
                maxCond = 10^(max(obj.N/5, 1));
                if nargin < 2
                    K = 1e2; % default: 100 samples
                end
            end
            
            % sample K combinations of 2N-3 arcs
            allCombs = obj.N * (obj.N-1) / 2;
            selCombs = 2 * obj.N - 3;
            randCombination = zeros(selCombs, K);

            % iterate over all sampled combinations
            doaCombinations = nchoosek(1:obj.N, 2);
            S_barOpt = zeros(selCombs, selCombs, K);
            [conditioning, ~] = deal(NaN(K, 1));
            for ki = 1:K
                
                % draw random sample from possible arc combinations
                randCombination(:, ki) = sort(randsample(allCombs, selCombs));
                arcSelection = doaCombinations(randCombination(:, ki), :);
                
                % build covariance matrix for this selection
                S_barOpt(:, :, ki) = obj.getCovarianceMatrix(arcSelection);
                
                % get conditioning of covariance
%                 [U, S, V] = svd(S_barOpt(:, :, ki));
                conditioning(ki) = cond(S_barOpt(:, :, ki));
                
%                 % get Mahalanobis distance
%                 if conditioning(ki) < 1e8
%                     phi_bar = obj.getPhi_bar(arcSelection);
%                     mahal(ki) = phi_bar' * (V/S*U') * phi_bar;
%                 end
                
                % exit condition: conditioning below 100
                if conditioning(ki) < maxCond
                    break
                end
            end
            
            [minCond, KI] = min(conditioning);
            
            arcs = doaCombinations(randCombination(:, KI), :);
            S_bar = S_barOpt(:, :, KI);
        
        end
        
        
        function S_bar_sel = getCovarianceMatrix(obj, gcaSelection)
            %obj.getCovarianceMatrix(gcaSelection) construct Covariance
            %matrix
            %   For a selection of arcs constructs the respective
            %   covariance matrix.
            
            % preallocate correlation matrix
            correlation = zeros(2*obj.N-3, 2*obj.N-3);
            
            % loop over DoAs
            for n = 1:obj.N
                % find arcs that use this DoA
                arcs = find(any(gcaSelection==n, 2));
                
                if numel(arcs) > 1
                    % get all combinations of arcs with this angle
                    gs = nchoosek(arcs, 2);
                    % gedefinet adjacent angles
                    toprow = gcaSelection(gs(:, 1), :)';
                    toprow = toprow(toprow ~= n);
                    botrow = gcaSelection(gs(:, 2), :)';
                    botrow = botrow(botrow ~= n);
                    % compute spherical cosine
                    sphCos = obj.sphericalCosine( ...
                        [toprow'; repmat(n, size(toprow')); botrow']);
                    % assign values in correlation matrix
                    correlation(sub2ind(size(correlation), gs(:, 1), gs(:, 2))) = ...
                        sphCos * obj.sigmas(n)^2;
                end
            end
            
            %  scale correlation with angles that are close to each other
            centralAngles = obj.gcas(sub2ind(size(obj.gcas), ...
                gcaSelection(:, 1), gcaSelection(:, 2)));
            normAnglesSq = centralAngles.^2 ./ sum(obj.sigmas(gcaSelection).^2, 2);
            scaleMatrix = diag(1-exp(-0.5*normAnglesSq));
            
            correlation = scaleMatrix * correlation * scaleMatrix;
            
            % construct covariance matrix
            S_bar_sel = diag(sum(obj.sigmas(gcaSelection).^2, 2)) ...
                + correlation + correlation';
            
        end
        
        
        function gcas = getGCAs(obj, gcaSelection)
            %obj.getGCAs(gcaSelection) computes vector of arcs
            %   For a specific selection of arcs, creates the vector of
            %   these arcs.
            %   
            %   Input:
            %   gcaSelection    M x 2 matrix of tuples of indices of 
            %                   begin/end doas for M arcs
            
            % select gcas from precomputed matrix
            gcas = obj.gcas(sub2ind(size(obj.gcas), ...
                gcaSelection(:, 1), gcaSelection(:, 2)));
        end
        
        
        function cosTheta = sphericalCosine(obj, doaIndices)
            %obj.sphericalCosine(doaIndices) Computes spherical cosine in
            %triangle defined by 3 doaIndices at middle index.
            %   
            %   Inputs:
            %   doaIndices  3 element vector of indices of doas
            %   
            %   Outputs:
            %   cosTheta    Cosine of spherical angle enclosed by the
            %               triangle
            
            
            % get involved great circle arcs
            gcaI = obj.gcas(sub2ind(size(obj.gcas), ...
                doaIndices, [0 1 0; 0 0 1; 1 0 0]*doaIndices));
            
            % compute spherical cosine
            cG = cos(gcaI);
            sG = sin(gcaI(1:2, :));
            cosTheta = (cG(3, :) - cG(1, :).*cG(2, :)) ./ sG(1, :) ./ sG(2, :);
        end
        
        function gca = computeGCAs(obj, gcaSelection, unitVectors)
            %obj.computeGCAs(gcaSelection, unitVectors) computes arcs
            %   For a given selection of doa combinations, computes great
            %   circle arcs between said arcs.
            %   If no matrix of arcs is passed uses obj.doaUnitVectors.
            %   
            %   Input:
            %   gcaSelection    M x 2 matrix of tuples of indices of 
            %                   begin/end doas for M arcs
            %   unitVectors     (optional) 3 x K matrix of unit vectors
            
            if nargin < 3
                unitVectors = obj.doaUnitVectors;
            end
            % compute gcas
            gca = squeeze(acos(dot(unitVectors(:, gcaSelection(:, 1), :), ...
                                   unitVectors(:, gcaSelection(:, 2), :))));
        end
        
        
    end
end

