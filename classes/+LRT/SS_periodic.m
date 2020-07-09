classdef SS_periodic < LRT.SS
    %class of simple vs. simple Likelihood Ratio of periodic measurements 
    %   Calculates a simple vs. simple likelihood ratio test for a set of 
    %   measured values, actual values, measurement standard deviation, etc.
    %   
    %   Implementation using dimensionality reduction to run hypotheses
    %   tests using N-1 measurements of the differences between satellite
    %   azimuths.
    %   Contains methods to exclude outliers to H0, outliers to H1.
    
    
    properties
        p               % period of measurement
    end
    
    properties (Dependent)
        phi_bar         % adjusted truth value under H0
        y_bar           % adjusted measurement
        err_bar         % adjusted measurement error
        R_bar           % adjusted covariance matrix
    end
    
    
    methods
        function obj = SS_periodic(period, varargin)
            %SS_periodic(period, mu0, mu1, sigma, measured, SV)
            %   Constructs an instance of the class for specific period,
            %   expected values, measurements, measurement standard
            %   deviations and satellite names.
            
            % call superclass constructor
            obj = obj@LRT.SS(varargin{:});
            
            if nargin > 0
                obj.p = period;
            else
                obj.p = NaN;
            end
            
        end
        
        % ------------ GET METHODS ---------------------
        function phi_bar = get.phi_bar(obj)
            %phi_bar Adjusted measurement truth
            %   Calculates the adjusted measuremend truth under
            %   dimensionality reduction.
            
            s2u = obj.sats2use; % grab once for speed
            
            % perform dimensionality reduction
            phi_bar = obj.wrap(obj.reductionMatrix * obj.mu0(s2u));
            if ~isempty(phi_bar)
                % choose optimal wrapping
                phi_bar = obj.minMahalanobis(phi_bar);
                % run twice for robustness / two steps might still reduce
                % it
                phi_bar = obj.minMahalanobis(phi_bar); 
            end
        end
        
        function y_bar = get.y_bar(obj)
            %obj.y_bar Adjusted measurement
            %   Calculates the adjusted measuremend  under
            %   dimensionality reduction.
            
            s2u = obj.sats2use; % grab once for speed
            
            % perform dimensionality reduction
            A = obj.reductionMatrix(obj.mu0(s2u));
            
            % wrap to [-obj.p/2 obj.p/2]
            y_bar = obj.wrap(A * obj.y(s2u));
            
%             if ~isempty(y_bar)
%                 y_bar = obj.minMahalanobis(y_bar);
%             end
            
        end
        
        function err_bar = get.err_bar(obj)
            %obj.DDerr
            %   Measured-expected double differences after division by 10 
            %   in units of cycles. Represents the double differences
            %   error after removing the integer ambiguity.
            
            err_bar = obj.H0error(obj.y);
        end
        
        function Rb = get.R_bar(obj)
            %obj.R_bar Covariance matrix of DD measurement
            %   Computes the covariance matrix of the double difference
            %   measurements in units of cycles^2.
            
            % get reduction matrix
            A = obj.reductionMatrix;
            
            % get DD covariance matrix in units of cycles^2
            Rb = A * diag(obj.sigma(obj.sats2use).^2) * A';
        end
        
        % ------------ HELPER METHODS -------------------
        function v = wrap(obj, values)
            %obj.wrap(values) Wraps values to [-0.5 0.5]*obj.p
            %   Wraps values around 0 given the periodicity p. Operates on
            %   an element by element basis.
            v = (values/obj.p - round(values/obj.p)) * obj.p;
        end
        
        function H0err = H0error(obj, y)
            %H0error(y) calculates error of measurement under H0
            %   H0err = obj.H0error(y)
            %   Calculates the distance of a measurement from mu0 under the
            %   dimensionality reduction scheme.
            
            s2u = obj.sats2use; % grab once for speed
            
            % check passed vector size
            if size(y, 1) > sum(s2u)
                y = y(s2u, :);
            end
            
            if size(y, 1) < 2 % need min 2 measurements
                H0err = [];
            else
                % calculate distances from mu_0, mu_1
                H0err = obj.minMahalanobis(obj.minMahalanobis( ...
                    obj.wrap(obj.reductionMatrix * (y - obj.mu0(s2u)))));
            end
        end
        
        function H1err = H1error(obj, y)
            %H1error(y) calculates error of measurement under H1
            %   H1err = obj.H1error(y)
            %   Calculates the distance of a measurement from mu1 under the
            %   dimensionality reduction scheme.
            
            s2u = obj.sats2use; % grab once for speed
            
            % check passed vector size
            if size(y, 1) > sum(s2u)
                y = y(s2u, :);
            end
            
            if size(y, 1) < 2
                H1err = [];
            else
                % calculate distances from mu_0, mu_1
                H1err = obj.minMahalanobis(obj.minMahalanobis( ...
                    obj.wrap(obj.reductionMatrix * (y - obj.mu1(s2u)))));
            end
        end
        
        function p_yH0 = getP_yH0(obj, y)
            %Calculates the conditional probability of y given H0
            %   p_yH0 = obj.getP_yH0(y)
            %   Calculates the conditional probability of the measured
            %   azimuths given the nominal Hypothesis by evaluating the
            %   multivariate normal distribution of proper covariance 
            %   at the differences in measurement errors.
            %   
                        
            if sum(obj.sats2use) > 1
                p_yH0 = mvnpdf(obj.H0error(y)', [], obj.R_bar);
            else
                p_yH0 = NaN;
            end
        end
        
        function p_yH1 = getP_yH1(obj, y)
            %Calculates the conditional probability of y given H1
            %   p_yH1 = obj.getP_yH1(y)
            %   Calculates the conditional probability of the measured
            %   azimuths given the spoofed Hypothesis by evaluating the
            %   multivariate normal distribution of proper covariance 
            %   at the differences in measurement errors.
            %   
            
            nSats = sum(obj.sats2use);
            
            if nSats > 1
                % get Mahalanobis distances of a couple of options                
                [~, ~, mahalDists] = obj.minMahalanobis(obj.H1error(y));
                
                % average over a range of fits to different Gaussians
                likelihoods = (2*pi)^(-(obj.N-1)/2) ...
                    * det(obj.R_bar)^(-1/2) * exp(-1/2*mahalDists);
                
                % take average over all Gaussians that contribute
                % average over all if none contribute
                l2take = likelihoods > 10^(-nSats);
                p_yH1 = sum(likelihoods, 2) ...
                    ./ (sum(l2take, 2) + ~any(l2take, 2)*size(l2take, 2));
            else
                p_yH1 = NaN;
            end
        end
        
        function lLy = getLogLambda(obj, y)
            %obj.getLogLambda(y) Computes the log Lambda of a measurement y
            %   For a given measurement y, calculates the resulting value
            %   of the decision variable log Lambda. If no argument is
            %   passed, returns the saved value of obj.logLambda.
            
            if nargin < 1
                y = obj.y;
            end
            
            if sum(obj.sats2use) < 2
                lLy = NaN;
            else
                RIhalf = chol(eye(sum(obj.sats2use)-1) / obj.R_bar);

                lLy = - sum((RIhalf * obj.H0error(y)).^2, 1) / 2 ...
                      + sum((RIhalf * obj.H1error(y)).^2, 1) / 2;
                % potentially switch this to log(p_yH0) - log(p_yH1) for
                % averaged p_yH1 computation
            end
        end
        
        function A = reductionMatrix(obj, angles)
            %reductionMatrix(obj, angles) dimensionality reduction matrix
            %   The difference between angles is used in the hypothesis
            %   tests. This matrix calculates these differences and reduces
            %   the dimension by 1 for the selected satellites. We build
            %   the matrix for angles sorted according to increasing
            %   azimuth.
            if nargin < 2
                angles = obj.mu0(obj.sats2use);
            end
            if any(angles)
                I = length(angles) - 1;
                % sort for increasing azimuths for smoother behavior
                [~, sI] = sort(angles);
                
                A = full(sparse([1:I, 1:I], [sI(1:end-1) sI(2:end)], ...
                        [-ones(1, I), ones(1, I)], I, I+1));
                
                % alternatively banded diagonal matrix
%                 B = full(spdiags([ones(I, 1), -ones(I, 1)], [0 1], I, I+1));
            else
                A = [];
            end
        end
        
        function [mu, Sigma] = nominalDistribution(obj)
            %obj.nominalDistribution calculates mean and covariance
            %   of the log-distribution expected under H_0.
            %   Overwrites same method of superclass due to dimensionality
            %   reduction.
            
            phi_b = obj.phi_bar;
            
            if ~isempty(phi_b)
                mu = 0.5 * phi_b' / obj.R_bar * phi_b; % mean of Lambda
            else
                mu = NaN;
            end
            
            if nargout > 1
                Sigma = 2 * mu; % Variance of Lambda
            end
        end
        
        function SVstrings = SVdifferences(obj)
            %SVstrings = SVdifferences(obj)
            %   Creates cell array of strings indicating dimensionality
            %   reduction differences taken.
            
            SVs = obj.SVnames(obj.sats2use);
            A = obj.reductionMatrix;
            
            S = sum(obj.sats2use) - 1;
            SVstrings = cell(S, 1);
            
            for s = 1 : S
                SVstrings{s} = strjoin({SVs{A(s, :) == 1}, ...
                    '-', SVs{A(s, :) == -1}});
            end
            
        end
        
        function [optValues, minDist, allDist] = minMahalanobis(obj, values, constr)
            %[values, minDist] = obj.minMahalanobis(values)
            %   Wrap values for lowest mahalanobis distance for given 
            %   covariance matrix R_bar and wrapping limit.
            %   Works on matrix of multiple values, expecting column
            %   vectors.
            %   Includes the option to set a minimum value as constraint
            %   for the Mahalanobis distance.
            
            if nargin < 3
                constr = 0;
            end
            
            [I, K] = size(values);
            
            RinvC = chol(eye(I) / obj.R_bar); % precompute for speed
            
            % normalize by period
            valuesN = values / obj.p;
            
            % build alternative options
            
            haveSigns = sign(values);
            
            % deeper sign offsets
            diagMat = diag(ones(I, 1)) + diag(ones(I-1, 1), 1);
            diagMat(end, 1) = 1;
            
            % collect all values in matrix to choose best
            valuesDiags = zeros(I, K*2*I);
            for i = 1:I
                valuesDiags(:, (i-1)*K+1 : i*K) = ...
                    - [zeros(i-1, K); haveSigns(i, :); zeros(I-i, K)];
                valuesDiags(:, (i+I-1)*K+1 : (i+I)*K) = ...
                    - haveSigns .* diagMat(:, i);
            end
                       
            % check for optimal distribution of alternating signs
            wantSigns = (-1).^(0:I-1)' * haveSigns(1, :);
            
            abovePover4 = abs(values) > obj.p/4;
            posOffset = abovePover4 .* (wantSigns - haveSigns);
            negOffset = abovePover4 .* (- wantSigns - haveSigns);
            
            % concatenate options in I x K(2*I+3) matrix
            allValuesN = [zeros(I, K), valuesDiags, ...
                          posOffset, negOffset] + repmat(valuesN, 1, 2*I+3);
            
            % finally calculate mahalanobis distances, find min distance
            mahalSums = sum((RinvC * allValuesN).^2, 1);
            allMahals = reshape(mahalSums, K, 2*I+3);

            [minD, Imi] = min(allMahals, [], 2);
            
            % select optimal values, unnormalize
            optValues = allValuesN(:, (1:K) + (Imi'-1)*K) * obj.p;
            
            % unnormalize for mahalanobis distance
            if nargout > 1
                minDist = minD * obj.p^2;
                if nargout > 2
                    allDist = allMahals * obj.p^2;
                end
            end
            
        end
        
        % -------- PLOTTING METHODS -------------------
        
        function f = plotAdjMeas(obj, satellites)
            %f = obj.plotAdjMeas(satellites)
            %   Plots the adjusted measured values and their uncertainty as
            %   errorbar plot against the true values for a selection of
            %   satellites. Defaults to obj.sats2use if no satellite
            %   selection is specified.
            
            if nargin < 2
                satellites = obj.sats2use;
            end
            
            A = obj.reductionMatrix(obj.mu0(satellites));
            
            
            % check if hold on is active
            HoldFlag = ishold;
            
            % plot errorbar plot of measurements
            errorbar(obj.y_bar / obj.p, ...
                     sqrt(diag(obj.R_bar)) / obj.p, ...
                     'b+')
            f = gca;
            hold on; grid on;
            % plor actual azimuths
            plot(obj.phi_bar / obj.p, 'r.', 'MarkerSize', 25)
            
            % create x-tick labels
            
            SVs = find(satellites);
            S = sum(satellites) - 1;
            diffStrings = cell(S, 1);
            for s = 1 : S
                diffStrings{s} = strjoin({num2str(SVs(A(s, :) == 1)), ...
                    '-', num2str(SVs(A(s, :) == -1))});
            end
            f.XTick = 1:S;
            f.XTickLabel = diffStrings;
            
            % label
            xlabel('Satellites', 'FontSize', 16, 'Interpreter', 'latex')
            ylabel('Fractions of period', ...
                'FontSize', 16, 'Interpreter', 'latex')
            legend('Measured, $2\sigma$', 'Actual', ...
                'Location', 'best', 'FontSize', 16, 'Interpreter', 'latex')
            
            
            if ~HoldFlag
                hold off;
            end
        end
    end
end

