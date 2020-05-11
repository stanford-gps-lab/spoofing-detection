function r = qMult(p, q)
%qMult(p, q)
%   Quaternion multiplication p * q.
%   Follows the convention that multiplication is right to left.
%   
%   Expects p and q to have the same dimensions. They can be 4 x N or N x 4
%   arrays.
%   
%   Input:
%   p   array<double> of quaternions. Can be 4 x N or N x 4.
%   q   array<double> of quaternions. Must be same dimensions as p.
%   
%   Output:
%   r   array<double> of quaternions. Same dimensions as p.
%   
%   Based on the formulas in "Representing Attitude: Euler Angles, Unit 
%   Quaternions and Rotation Vectors" by James Diebel, 20 Oct. 2006.
%   
%   Written by Fabian Rothmaier, 2019

% ensure working with column vectors of quaternions
if size(q, 1) ~= 4
    flipDim = true;
    q = q';
    p = p';
else
    flipDim = false;
end

% compute multiplication for each quaternion pair
r1 = sum([p(1, :); -p(2:4, :)] .* q, 1);
r2 = sum([p(2, :); p(1, :); p(4, :); -p(3, :)] .* q, 1);
r3 = sum([p(3, :); -p(4, :); p(1, :); p(2, :)] .* q, 1);
r4 = sum([p(4, :); p(3, :); -p(2, :); p(1, :)] .* q, 1);

r = [r1; r2; r3; r4];

% revert dimensions back if necessary
if flipDim
    r = r';
end

end

