function sigmaDelta = deltaSigmaVal(CN0, noiseModel)
    %Compute the standard deviation of the delta metric on the tracking tap
    %as a function of C/N0. Is set up to represent a Novatel GIII receiver.

    if nargin < 2
        noiseModel = 'simple';
    end

    if strcmp(noiseModel, 'Betz')

        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %Use model from Betz and Kolodziejski, 2000

        Tc = 1/GPSconstants.chippingRate('L1'); % chip length in s
        T = 0.02; % pre-detection integration time 20 ms

        B_L = 15; %4;

        D = 0.1032; % spacing between tracking taps in chips
        b = 24e6 * Tc; % For 24 MHz Bandwidth
        Db = D * b;

        if Db >= pi
            SigmaNorm = B_L * (1-0.5*B_L*T) ./ 10.^(CN0/10) / 2 * D;

        elseif Db > 1
            SigmaNorm = B_L * (1-0.5*B_L*T) ./ 10.^(CN0/10) / 2 ...
                * (1/b + b/(pi-1)*(D-1/b)^2);
        else % Db <= 1
            SigmaNorm = B_L * (1-0.5*B_L*T) ./ 10.^(CN0/10) / 2 / b;
        end

        sigmaDelta = sqrt(2) * sqrt(SigmaNorm);

    elseif strcmp(noiseModel, 'simple')
        % use simple model
        sigmaDelta = 0.1 ./ (CN0 - 30);
    else
        error('Invalid noise model.');
    end

end