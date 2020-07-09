classdef acfDeltaMetric < LRT.SC_GLRT
    %acfDeltaMetric Likelihood ratio test of acf Delta metric
    %   Performs a simple vs. composite hypothesis test in form of a
    %   Generalized Likelihood Ratio Test (GLRT) for delta metrics of
    %   correlator tabs of the autocorrelation function (acf).
    %   
    %   Calls its superclass LRT.SC_GLRT when constructing the object.
    %   
    %   @params:
    %   
    
    properties
        tabSigma    % Covariance matrix of correlator tabs
        Sigma       % Covariance matrix of delta metrics
    end
    
    methods
        
        function obj = acfDeltaMetric(sigmas, spacing, varargin)
            %acfDeltaMetric(sigmas, spacing, varargin)
            %   Constructs a GLRT object for N acf delta metrics.
            %   
            %   @params:
            %   sigmas  <array> 2Nx1 of stdev. values of each correlator
            %   spacing <array> 2Nx1 of spacing in chips of each correlator 
            %           from the prompt
            %   
            
            % check inputs, fill empty inputs
            if nargin < 2
                if nargin < 1
                    sigmas = [];
                end
                spacing = ones(size(sigmas));
            end
            
            if length(sigmas) ~= length(spacing)
                error('sigmas and spacing vector must have equal length.')
            end
            if mod(length(sigmas), 2) ~= 0
                error('Delta metric needs even number of correlator tabs.')
            end
            
            % calculate Cov matrix of correlator tabs
            tabSigma = (sigmas(:) * sigmas(:)') ...
                .* LRT.acfDeltaMetric.Rcorr(spacing - spacing');
            
            % calculate Cov matrix of delta metrics
            N = length(sigmas) / 2;
            A = [eye(N), -fliplr(eye(N))];
            Sigma = A*tabSigma*A';
            SigmaInv = eye(N) / Sigma;
            
            
            % call superclass constructor
            obj@LRT.SC_GLRT(SigmaInv, zeros(size(sigmas)), 0, varargin{:});
            
            % set remaining parameters
            obj.tabSigma = tabSigma;
            obj.Sigma = Sigma;
            
        end
        
        
    end
    
    methods (Static = true, Sealed = true)
        function R = Rcorr(tau)
            %y_b = obj.prBias(xBias, satellites)
            %   Calculates the pseudorange bias to achieve a certain state
            %   bias. Can be called with a logical vector indicating which
            %   pseudoranges are modified.
            
            R = max(1 - abs(tau), 0);
            
        end
        
    end
end

