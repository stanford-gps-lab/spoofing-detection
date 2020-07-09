classdef SC_GLRT < LRT.GLRT
    %class of Generalized Likelihood Ratio Test for well simple H0 and
    %composite H1
    %   Calculates a likelihood ratio of a Gaussian variable. Hypothesis H0
    %   is defined by known mean and covariance, hypothesis H1 is
    %   defined by known covariance and unknown mean.
    %   
    
    properties
        mu0             % expected value of measurement under H0
        p               % Number of dof less than number of measurements
    end
    
    properties (Dependent)
        subsetCount     % number of considered subsets during iteration
        dof             % degrees of freedom / number of used satellites
    end
    
    
    methods
        function obj = SC_GLRT(SigmaInv, mu0, p, varargin)
            %SC_GLRT(SigmaInv, mu0, p, y, SVnames)
            %   Constructs an instance of the class.
            %   
            %   Inputs:
            %       SigmaInv    <matrix> inverse of measurement covariance 
            %                   matrix. Defaults to empty matrix.
            %       mu0         [optional]<vector> measurement mean under 
            %                   H0. Defaults to zeros(size(SigmaInv, 1), 1)
            %       p           [optional]<int> number of deg. of freedom 
            %                   less than number of measurements. Defaults 
            %                   to 0.
            %       y           [optional]<vector> measurements
            %       SVnames     [optional]<cell><string> SV names
            
            if nargin == 0
                SigmaInv = [];
            end
            
            % call superclass constructor
            obj@LRT.GLRT(SigmaInv, varargin{:});
            
            % set remaining properties
            if nargin < 3
                p = 0;
                if nargin < 2
                    mu0 = zeros(size(SigmaInv, 1), 1);
                end
            end
            
            obj.mu0 = mu0(:);
            obj.p = p;
            
        end
        
        % ------------ GET METHODS ---------------------
        
        function k = get.dof(obj)
            %k = obj.dof Compute degrees of freedom
            %   Computed degrees of freedom = number of used measurements
            %   minus number of states
            k = sum(obj.sats2use) - obj.p;
        end
        
        function ssc = get.subsetCount(obj)
            %ssc = obj.subsetCount
            %   Computed the number of considered satellite subsets to be
            %   used to inflate the maximum false alert probability
            %   constraint.
            
%             % conservative approach: sum of binomial coefficients
%             k = sum(obj.consideredSats):obj.N;
%             ssc = sum( ...
%                 factorial(obj.N) ./ factorial(k) ./ factorial(obj.N-k));
            
            % approach of counting only considered sets
            K = sum(obj.consideredSats);
            ssc = 1 + 0.5 * (obj.N^2 + obj.N - (K^2 + K));
            
        end
        
        % ------------ HELPER METHODS -------------------
        
        function H0err = H0error(obj, y)
            %obj.H0error(y) calculate the difference between a measurement y
            %and its mean under nominal conditions.
            %   
            if size(y, 1) < length(obj.mu0)
                H0err = y - obj.mu0(obj.sats2use);
            else
                H0err = y - obj.mu0;
            end
        end
        
        function H1err = H1error(obj, y)
            %obj.H1error(y) calculate the difference between a measurement y
            %and its mean under spoofed conditions.
            %   Resulting error is always zero, the mean under H1 is equal
            %   to the measurement as it is the MLE mean.
            
            H1err = zeros(obj.N, size(y, 2));
        end
        
        function p_H0 = getP_logLambdaH0(obj, logLambda)
            %p_H0 = obj.getP_logLambdaH0(logLambda) evaluate pdf | H0
            %   Evaluates the pdf of log Lambda for passed values of
            %   logLambda.
            
            p_H0 = 2 * chi2pdf(-2*logLambda, obj.dof);
            
        end
        
        function p_H1 = getP_logLambdaH1(obj, logLambda, delta)
            %p_H1 = obj.getP_logLambdaH1(logLambda, delta) evaluate pdf | H1
            %   Evaluates the pdf of log Lambda for passed values of
            %   logLambda and noncentrality parameter delta.
            
            if nargin < 3
                delta = NaN;
            end
            p_H1 = 2 * ncx2pdf(-2*logLambda, obj.dof, delta);
            
        end
        
        function obj = excludeMaxErr(obj)
            %obj.excludeMaxVal Determine and exclude the measurement with
            %the largest error.
            %   Finds and excludes the satellite with the largest
            %   measurement error. Exclusion is shown by setting the
            %   respective value in obj.excludedSats = true.
            
            % find largest error
            [~, maxSat] = max(obj.y(obj.sats2use) - obj.mu0(obj.sats2use));
            % get indices of used satellites
            usedIndices = find(obj.sats2use);
            
            % exclude outlier
            obj.excludedSats(usedIndices(maxSat)) = true;
        end
        
        function gamma = threshold(obj, P_FAmax, dof)
            %threshold(P_FAmax, dof) calculates the alarm threshold
            %   gamma = obj.threshold(P_FAmax, dof)
            %   Calculates the Executive Monitor (EM) alarm threshold that
            %   satisfies the specified false alert probability. Can be
            %   called with a specific number of degrees of freedom.
            
            if nargin < 3
                dof = obj.dof;
            end
            % calculate threshold
            gamma = - 1/2 * chi2inv(1-P_FAmax, dof);   
        end
        
        function lambda = chi2statistic(obj, y_r)
            %lambda = obj.chi2statistic(y_r)
            %   Calcualte the chi2 statistic of the measurement y_r.
            
            lambda = sum((obj.Phalf * obj.H0error(y_r)).^2, 1);
        end
        
        function p = power(obj, P_FAmax, lambda)
            %power(P_FAmax, lambda) calculates the GLRT test power
            %   p = obj.power(P_FAmax, lambda)
            %   Calcualtes the power of the GLRT to detect an attack of
            %   signals for a noncentrality parameter under H1, lambda.
            
            % calculate power
            p = ncx2cdf(-2*obj.threshold(P_FAmax), obj.dof, ...
                        lambda, 'upper');   
        end
        
        % ------------ PLOTTING METHODS -------------------
        
        function f = plotLambdaDistribution(obj, varargin)
            %f = obj.plotLambdaDistribution(P_FAmax)
            %   Plot the Lambda decision space with nominal and spoofed
            %   distribtions, thresholds and measured slot.
            
            % use superclass method
            f = plotLambdaDistribution@LRT.GLRT(obj, varargin{:});
            xlim([-inf, 0]);
        end
    end
    
    
end

