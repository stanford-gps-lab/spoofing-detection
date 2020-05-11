function [azMeas, elMeas] = doa2Dnoise(azTrue, elTrue, sigma)
%[azMeas, elMeas] = doa2Dnoise(azTrue, elTrue, sigma)
%   Adds noise to 2D Direction of Arrival (DoA) measurements. Noise is
%   added in the form of a rotation in an arbitrary direction of the vector
%   by a Normally distributed magnitude.
%   
%   @params:
%   azTrue  <double>matrix of size N x K true azimuth values
%   elTrue  <double>matrix of size N x K true elevation values
%   sigma   <double>array of standard deviation of measurements
%           Can be size N x 1 or N x K.
%   
%   @out:
%   azMeas  <double>matrix of size N x K azimuth values of DoA measurements
%   elMeas  <double>matrix of size N x K elevation of DoA measurements

[N, K] = size(azTrue);

if size(sigma, 1) ~= N
    error('First dimensions of azimuth, sigma need to be of equal length.')
end

[elMeas, azMeas] = deal(zeros(N, K));
for n = 1:N
    % convert true direction to quaternions
    qTrue = euler2q([zeros(1, K); elTrue(n, :); azTrue(n, :)]);
    
    % if only spoofed measurements are generated
%     qTrue = euler2q([zeros(1, K); repmat(elSpoof, 1, K); repmat(azSpoof, 1, K)]);
    
    % generate randomized noise quaternion
    delta = sigma(n, :) .* randn(1, K); % magnitude of error
    alpha = 2*pi * rand(1, K); % direction of error
    qNoise = [cos(delta/2); ...
              sin(delta/2) .* [zeros(1, K); sin(alpha); cos(alpha)]];
    
    % compute measurement quaternion
    qMeas = qMult(qNoise, qTrue);
    
    % convert back to euler angles
    eulerMeas = q2euler(qMeas);
    
    % extract az, el
    elMeas(n, :) = eulerMeas(2, :);
    azMeas(n, :) = mod(eulerMeas(3, :), 2*pi);
end

% adjust for +90 el jump
azMeas(elMeas > pi/2) = mod(azMeas(elMeas > pi/2) + pi, 2*pi);
elMeas(elMeas > pi/2) = pi - elMeas(elMeas > pi/2);


end

