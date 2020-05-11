function q = euler2q(eulerAngles)
%euler2q(eulerAngles)
%   Calculates quaternions q from Euler angles roll, pitch, yaw. Returns
%   angles in a vector of dimensions similar to dimensions of eulerAngles.
%   Operates on each column of eulerAngles, unless the number of rows of 
%   eulerAngles is not 3.
%   Input:
%   eulerAngles array<double> of calculated Euler Angles. Can be 3 x N or 
%               N x 3.
%   
%   Output:
%   q           array<double> of quaternions. In dimensions
%               similar to q. That is, if q is a 4 x N matrix, then
%               eulerAngles is a 3 x N matrix and vice versa.
%   
%   Based on the formulas in "Representing Attitude: Euler Angles, Unit 
%   Quaternions and Rotation Vectors" by James Diebel, 20 Oct. 2006.
%   
%   Written by Fabian Rothmaier, 2019

% ensure working with column vectors of Euler angles
if size(eulerAngles, 1) ~= 3
    flipDim = true;
    seA = sin(eulerAngles' / 2);
    ceA = cos(eulerAngles' / 2);
else
    flipDim = false;
    seA = sin(eulerAngles / 2);
    ceA = cos(eulerAngles / 2);
end
q = [prod(ceA, 1) + prod(seA, 1); ...
     -ceA(1, :) .* seA(2, :) .* seA(3, :) + ...
        seA(1, :) .* ceA(2, :) .* ceA(3, :); ...
     ceA(1, :) .* seA(2, :) .* ceA(3, :) + ...
        seA(1, :) .* ceA(2, :) .* seA(3, :); ...
     ceA(1, :) .* ceA(2, :) .* seA(3, :) - ...
        seA(1, :) .* seA(2, :) .* ceA(3, :)];

% transpose again to original dimensions if necessary
if flipDim
    q = q';
end

end

