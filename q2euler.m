function eulerAngles = q2euler(q)
%q2euler(q)
%   Calculates Euler angles roll, pitch, yaw from quaternions q. Returns
%   angles in a vector of dimensions similar to dimensions of q.
%   Operates on each column of q, unless the number of rows of q is not 4.
%   
%   Input:
%   q           array<double> of quaternions. Can be 4 x N or N x 4.
%   
%   Output:
%   eulerAngles array<double> of calculated Euler Angles. In dimensions
%               similar to q. That is, if q is a 4 x N matrix, then
%               eulerAngles is a 3 x N matrix and vice versa.
%   
%   Based on the formulas in "Representing Attitude: Euler Angles, Unit 
%   Quaternions and Rotation Vectors" by James Diebel, 20 Oct. 2006.
%   
%   Written by Fabian Rothmaier, 2019

% ensure working with column vectors of quaternions
if size(q, 1) ~= 4
    flipDim = true;
    q = q';
else
    flipDim = false;
end
eulerAngles = [atan2(2.*q(3, :).*q(4, :) + 2.*q(1, :).*q(2, :), ...
                     q(1, :).^2 + q(4, :).^2 - q(2, :).^2 - q(3, :).^2); ...
               -asin(2.*q(2, :).*q(4, :) - 2.*q(1, :).*q(3, :)); ...
               atan2(2.*q(2, :).*q(3, :) + 2.*q(1, :).*q(4, :), ...
                     q(1, :).^2 + q(2, :).^2 - q(3, :).^2 - q(4, :).^2)];

% transpose again to original dimensions if necessary
if flipDim
    eulerAngles = eulerAngles';
end

end

