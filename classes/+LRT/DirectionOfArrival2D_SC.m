classdef DirectionOfArrival2D_SC < LRT.SC_GLRT
    %DirectionOfArrival2D_SC
    %   Performs a simple vs. composite hypothesis test in form of a
    %   Generalized Likelihood Ratio Test (GLRT) for 2 dimensional
    %   Direction of Arrival (DoA) measurements.
    %   
    %   Calls its superclass LRT.SC_GLRT when constructing the object.
    %   
    %   @params:
    %   az      satellite azimuths
    %   el      satellite elevations
    %   ephemUV ephemeris DoA unit vectors in local coordinate frame
    
    properties
        az      % satellite azimuths (rad)
        el      % satellite elevation (rad)
        ephemUV % unit vectors in local coordinate frame of satellite DoAs
    end
    
    methods
        
        function obj = DirectionOfArrival2D_SC(az, el, weights, varargin)
            %pseudorangeResiduals(az, el, weights, y, SVnames)
            %   Constructs a GLRT object for pseudorange residuals.
            %   
            
            % generate weighting matrix
            W = diag(weights.^2);
            
            % call superclass constructor
            obj@LRT.SC_GLRT(W, [az; el], 3, varargin{:});
            
            % set additional properties
            obj.az = az;
            obj.el = el;
            obj.ephemUV = azel2enu(az, el);
        end
        
        function H0err = H0error(obj, y_r)
            %obj.H0error(y) calculate the difference between a measurement y
            %   and its mean under nominal conditions.
            %   Overwrites the superclass method to allow for attitude
            %   computation first.
            
            % compute mean first
            y_az = y_r(1:obj.N, :);
            y_el = y_r(obj.N+1:end, :);
            
            yUVs = azel2enu(y_az, y_el);
            
            % calculate H_0 error for every measurement
            H0err = zeros(obj.N, size(y_r, 2));
            for yi = 1:size(y_r, 2)
                [U, ~, V] = svd(yUVs(:, :, yi) * obj.ephemUV');
                R = U * diag([1, 1, det(U*V')]) * V';
                
                H0err(:, yi) = acos(diag(obj.ephemUV' * R' * yUVs(:, :, yi)));
            end
        end
        
    end
    
    methods (Static, Access = protected)
        
        function enu = azel2enu(azim, elev)
        %enu = azel2enu(azim, elev)
        %   Converts azimuth and elevation values to unit vectors in 
        %   east-north-up coordinate frame.
        %   Azimuth and elevation inputs must be vectors of equal size.
        %   
        %   Output matrix is 3 x size(azim)
        %   
        %   @params:
        %   azim    matrix of azimuth values in [rad]
        %   elev    matrix of elevation values in [rad]
        %   
        %   @out:
        %   enu     matrix of unit vectors
        %   

        % check input dimensions
        if size(azim) ~= size(elev)
            error('azel2enu must be called with vectors of equal size.')
        end

        % build unit vector elements
        enu1 = sin(azim) .* cos(elev);
        enu2 = cos(azim) .* cos(elev);
        enu3 = sin(elev);

        % construct result matrix
        enu = zeros([3, size(azim)]);
        enu(1, :) = reshape(enu1, 1, numel(enu1));
        enu(2, :) = reshape(enu2, 1, numel(enu2));
        enu(3, :) = reshape(enu3, 1, numel(enu3));

        end
        
    end
    
end

