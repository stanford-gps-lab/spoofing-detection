classdef pseudorangeResiduals < LRT.SC_GLRT
    %pseudorangeResiduals Likelihood ratio test of pseudorange residuals
    %   Performs a simple vs. composite hypothesis test in form of a
    %   Generalized Likelihood Ratio Test (GLRT) for pseudorange residuals.
    %   This is equivalent to a chi2 test of the residuals.
    %   
    %   Calls its superclass LRT.SC_GLRT when constructing the object.
    %   
    %   @params:
    %   az      satellite azimuths
    %   el      satellite elevations
    %   weights pseudorange weights (1/sigma)
    %   y       [optional] pseudorange updates
    %   SVnames [optional] satellite SV names
    
    properties
        az      % satellite azimuths (rad)
        el      % satellite elevation (rad)
        w       % pseudorange weights (1 / sigma)
        G       % geometry matrix
        W       % pseudorange weighting matrix (inverse covariance matrix)
        S       % position estimator (LS) x = S*y
    end
    
    methods
        
        function obj = pseudorangeResiduals(az, el, weights, varargin)
            %pseudorangeResiduals(az, el, weights, y, SVnames)
            %   Constructs a GLRT object for pseudorange residuals.
            %   
            
            % check inputs, fill empty inputs
            if nargin < 3
                if nargin < 2
                    if nargin < 1
                        az = [];
                    end
                    el = ones(size(az));
                end
                weights = ones(size(az));
            end
            
            % geometry matrix
            G = [-sin(az).*cos(el), ...
                 -cos(az).*cos(el), ...
                 -sin(el), ...
                 ones(size(az))];
            % weight matrix
            W = diag(weights.^2);
            
            % pseudorange residual covariance matrix
            WG = W*G;
            Sest = (G'*WG) \ WG'; % position estimator x = S*y
            P = W - WG*Sest;

            % call superclass constructor
            obj@LRT.SC_GLRT(P, zeros(size(az)), size(G, 2), varargin{:});
            
            % set additional properties
            obj.az  = az;
            obj.el  = el;
            obj.w   = weights;
            obj.G   = G;
            obj.W   = W;
            obj.S   = Sest;
            
        end
        
        function y_b = prBias(obj, xBias, satellites)
            %y_b = obj.prBias(xBias, satellites)
            %   Calculates the pseudorange bias to achieve a certain state
            %   bias. Can be called with a logical vector indicating which
            %   pseudoranges are modified.
            
            if nargin < 3
                satellites = true(size(obj.az));
            end
            
            % select estimator, weight matricies of satellite subset
            Sss = obj.S(:, satellites);
            Wss = obj.W(satellites, satellites);
            y_b = zeros(size(obj.az));
            
            % solve weighted least norm problem
            y_b(satellites) = Wss \ Sss' * ((Sss / Wss * Sss') \ xBias(:));
            
        end
        
    end
end

