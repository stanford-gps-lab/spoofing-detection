function enu = azel2enu(azim, elev)
%enu = azel2enu(azim, elev)
%   Converts azimuth and elevation values to unit vectors in east-north-up
%   coordinate frame.
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

