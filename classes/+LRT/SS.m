classdef SS < LRT.GLRT
    %class of simple vs. simple Likelihood Ratio
    %   Calculates a simple vs. simple likelihood ratio test for a set of 
    %   measured values, actual values, measurement standard deviation, etc.
    %   
    %   Contains methods to exclude outliers to H0, outliers to H1.
    
    properties
        mu0             % measurement truth | H0
        mu1             % measurement truth | H1
        sigma           % measurement standard deviation
    end
    
    properties (Dependent)
        normLogLambda   % normalized log of Lambda using nom. dist. method
    end
    
    methods
        function obj = SS(mu0, mu1, sigmas, varargin)
            %SS(mu0, mu1, sigmas, y, SVnames)
            %   Constructs an instance of the class for a mean under
            %   the nominal hypothesis, a mean under the alternate 
            %   hypothesis, measurements, measurement standard deviation
            %   and optionally satellite names.
            %   
            %   Supports measurement vectors. Treats each column as a
            %   separate measurement.
            
            if nargin < 2
                sigmas = [];
            end
            
            % call superclass constructor
            obj = obj@LRT.GLRT(diag(1./sigmas.^2), varargin{:});
            
            % set default values
            if nargin > 0
                obj.mu0 = mu0(:);
                obj.mu1 = mu1(:);
            else
                obj.mu0 = NaN;
                obj.mu1 = NaN;
            end
                        
            % read inputs
            if nargin > 2
                obj.sigma     = sigmas;
            end
            
            % update superclass definition of consideredSats
            obj.consideredSats  = obj.consideredSats & isfinite(obj.mu0);
            
        end
        % ------------ GET METHODS ---------------------
        
        
        function nLL = get.normLogLambda(obj)
            %normLogLambda The normalized log of Lambda
            %   The likelihood ratio Lambda has a certain expected
            %   normal distribution under nominal conditions. Normalizing
            %   creates a standard normally distributed variable that can
            %   easily be compared to the detection threshold.
            [mu, Sigma] = nominalDistribution(obj);
            nLL = (obj.logLambda - mu) ./ sqrt(Sigma);
        end
        
        % ------------ HELPER METHODS -------------------
        
        function H0err = H0error(obj, y)
            %obj.H0error(y) calculate the difference between a measurement 
            %y and its mean under nominal conditions.
            %   
            if size(y, 1) < length(obj.mu0)
                H0err = y - obj.mu0(obj.sats2use);
            else
                H0err = y - obj.mu0;
            end
        end
        
        function H1err = H1error(obj, y)
            %obj.H1error(y) calculate the difference between a measurement 
            %y and its mean under spoofed conditions.
            %   
            
            if size(y, 1) < length(obj.mu1)
                H1err = y - obj.mu1(obj.sats2use);
            else
                H1err = y - obj.mu1;
            end
        end
        
        function Rb = R(obj)
            %obj.R Covariance matrix of measurement
            %   Computes the covariance matrix of the measurements.
            
            Rb = diag(obj.sigma.^2);
        end
        
        function ssc = subsetCount(obj, K)
            %ssc = obj.subsetCount(K)
            %   Computed the number of considered satellite subsets to 
            %   adjust the maximum false alert probability. Counts all
            %   subsets down to the number of considered satellites K.
            %   Default for K is sum(obj.condiseredSats).
            
            if nargin == 1
                K = sum(obj.consideredSats, 1);
            end
            
            % conservative approach: sum of binomial coefficients
            k = K:obj.N;
            ssc = sum( ...
                factorial(obj.N) ./ factorial(k) ./ factorial(obj.N-k));
            
            % approach of counting only considered sets
%             ssc = 1 + 0.5 * (obj.N^2 + obj.N - (K^2 + K));
            
        end
        
        function [mu, Sigma] = nominalDistribution(obj)
            %obj.nominalDistribution calculates mean and covariance
            %   of the log-distribution expected under H_0.
            
            s2u = obj.sats2use;
            
            if any(s2u)
                % vectorized computation for speed
                RinvC = chol(eye(sum(s2u)) / obj.R(s2u, s2u));
                mu = 0.5 * sum((RinvC*(obj.mu0(s2u)-obj.mu1(s2u))).^2, 1);
            else
                mu = NaN;
            end
            
            if nargout > 1
                Sigma = 2 * mu; % Variance of Lambda
            end
        end
        
        function obj = detectMultipath(obj, k)
            %detectMultipath Determine the satellite affected by multipath
            %   Loops through all considered satellites and flags the
            %   satellites leading to the largest increase in the
            %   conditional probability of the nominal hypothesis. Only
            %   runs if no satellite has been flagged for multipath yet and
            %   if at least four satellites are considered.
            %   Flagging is done in the obj.multipath logical array.
            
            if nargin < 2
                k = 1;
            end
            
            % multipath exclusion loop:
            if sum(obj.sats2use) > 3
                nSats = length(obj.consideredSats(:, k));
                nominalLikelihood = NaN(nSats, 1);
                satellites2use = obj.sats2use(:, k); % precompute for speed
                for iSat = 1:nSats
                    if satellites2use(iSat) % if satellite used
                        % exclude satellite
                        obj.excludedSats(iSat, k) = true;
                        % save resulting conditional probability
                        nominalLikelihood(iSat) = obj.p_yH0;
                        % reset satellite
                        obj.excludedSats(iSat, k) = false;
                    end
                end
                % select highest probability, corresponding satellite
                [~, iMax] = max(nominalLikelihood);
                obj.excludedSats(iMax, k) = true;
            end
        end
        
        function [obj, iOutlier] = removeMaxOutlier(obj, k)
            %removeMaxOutlier(obj) removes largest outlier measurement
            %   iMax = removeMaxOutlier(obj, k)
            %   Removes the largest outlier satellite from being used among
            %   the kth measurement vector.
            %   The largest outlier is defined as the satellite who's
            %   removal leads to the largest increase in conditional
            %   probability of the spoofed hypothesis. This step resets the
            %   exclusion flag of the object.
            
            if nargin < 2
                k = 1;
            end
            
            % find satellite that makes the largest probability difference
            normLambdaPerSat = NaN(1, length(obj.consideredSats(:, k)));
            
            % calculate conditional probability for each removed satellite
            for iSat = 1:length(obj.consideredSats(:, k))
                if obj.consideredSats(iSat, k) % if satellite used
                    % reset multipath exclusion flag
                    obj.excludedSats(obj.excludedSats(:, k), k) = false;
                    % exclude satellite
                    obj.consideredSats(iSat, k) = false;
                    % reconsider multipath
                    obj = obj.detectMultipath(k);
                    % calculate respective normalized Lambda
                    normLambdaPerSat(iSat) = obj.normLogLambda;
                    % reset satellite
                    obj.consideredSats(iSat, k) = true;
                end
            end
            
            % pick satellite who's removal led to the smallest prob. ratio
            [~, iOutlier] = min(normLambdaPerSat);
            
            % remove selected satellite
            obj.consideredSats(iOutlier, k) = false;
            % reset exclusion flag
            obj.excludedSats(obj.excludedSats(:, k), k) = false;
            
        end
        
        function p_H0 = getP_logLambdaH0(obj, logLambda)
            %p_H0 = obj.getP_logLambdaH0(logLambda) evaluate pdf | H0
            %   Evaluates the pdf of log Lambda for passed values of
            %   logLambda.
            
            [mu, Sigma] = obj.nominalDistribution;
            
            p_H0 = normpdf(logLambda, mu, sqrt(Sigma));
            
        end
        
        function p_H1 = getP_logLambdaH1(obj, logLambda)
            %p_H1 = obj.getP_logLambdaH1(logLambda, delta) evaluate pdf | H1
            %   Evaluates the pdf of log Lambda for passed values of
            %   logLambda and noncentrality parameter delta.
            
            [mu, Sigma] = obj.nominalDistribution;
            
            p_H1 = normpdf(logLambda, -mu, sqrt(Sigma));
            
        end
        
        function ga = threshold(obj, P_FAmax)
            %threshold(P_FAmax) calculates the EM alarm threshold
            %   gamma = obj.threshold(P_FAmax)
            %   Calculates the Executive Monitor (EM) alarm threshold that
            %   satisfies the specified false alert probability.
            
            % get mean, variance of nominal distribution
            [mu, Sigma] = obj.nominalDistribution;
            % calculate threshold
            ga = norminv(P_FAmax, mu, sqrt(Sigma));   
        end
        
        function p = power(obj, P_FAmax)
            %power(P_FAmax) calculates the LR test power
            %   p = obj.power(P_FAmax)
            %   Calcualtes the power of the LR test to detect an attack of
            %   signals all coming from the same direction
            
            % get mean, variance of nominal distribution
            [~, Sigma] = obj.nominalDistribution;
            % calculate power
            x = sqrt(Sigma) + norminv(P_FAmax);
            if x > 5
                p = 1 - normcdf(x, 'upper');
            else
                p = normcdf(x);
            end
        end
        
        function pMD = missedDetectionP(obj, P_FAmax)
            %missedDetectionP(P_FAmax) calculates the LR P_MD
            %   p = obj.missedDetectionP(P_FAmax)
            %   Calcualtes the missed detection probability of the LR test 
            %   to detect an attack of signals all coming from the same 
            %   direction.
            
            % get mean, variance of nominal distribution
            [~, Sigma] = obj.nominalDistribution;
            % calculate quantile
            x = sqrt(Sigma) + norminv(P_FAmax);
            
            if x > 5
                pMD = normcdf(x, 'upper');
            else
                pMD = 1 - normcdf(x);
            end
        end
        
        % -------- PLOTTING METHODS -------------------
        
        
    end
end

