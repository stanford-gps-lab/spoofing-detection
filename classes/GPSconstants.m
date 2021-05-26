classdef GPSconstants < matlab.mixin.Copyable
    %GPSconstants class. Can be used as parent for a multitude of classes.
    %   Defines various constants and conversion factors.
    
    properties (Constant = true)
        % unit conversions
        m2ft = 1 / 0.3048;  % Meter to ft conversion
        NM2m = 1852;        % Nautical miles to meter conversion
        SM2m = 1609;        % Statue miles to meter conversion
        ms2kt = 36 / 18.52; % m/s to knots conversion
        d2r = pi/180;       % degree to radians
        r2d = 180/pi;       % radians to degree
        spw = 604800;       % seconds per week
        % common GPS constants
        c = 2.99792458e8;   % WGS-84 Speed of light in a vacuum (m/s)
        atomicStandard = 10.23e6;   % Hz atomic standard aboard satellite
    end
    
    methods
        
        function f = frequency(obj, freqName)
            %f = obj.frequency(frequencyName)
            % Calculates the center frequency in Hz.
            % Currently supports 'L1', 'L2' and 'L5'.
            f = obj.chippingRate(freqName) * obj.cyclesPerChip(freqName);
        end
        
        function c = chipLength(obj, freqName)
            %c = obj.chipLength(freqName)
            % Calculates the length of one chip in meter.
            c = obj.c / obj.chippingRate(freqName);
        end
        
        function l = lambda(obj, freqName)
            %l = obj.lambda(freqName)
            %   Calcualtes the wavelength lambda depending on the
            %   frequency.
            l = obj.c / obj.frequency(freqName);
        end
        
    end
    
    methods (Static = true)
                
        function cpc = cyclesPerChip(freqName)
            % Set cycles per code chip depending on frequency.
            % cyclesPerChip = 1540 for GPS L1 C/A
            %               = 1200 for GPS L2C
            %               = 115 for GPS L5
            
            if ~isa(freqName, 'char')
                % maybe it was passed as an integer?
                freqName = ['L', num2str(freqName)];
            end
            cycles = containers.Map({'L1', 'L2', 'L5'}, ...
                                    [1540, 1200, 115]);
            if cycles.isKey(freqName)
                cpc = cycles(freqName);
            else
                cpc = NaN;
            end
            
        end
        
        function cr = chippingRate(freqName)
            % chipping rate in chips per second
            
            if ~isa(freqName, 'char')
                % maybe it was passed as an integer?
                freqName = ['L', num2str(freqName)];
            end
            rates = containers.Map({'L1', 'L2', 'L5'}, ...
                                   [1, 1, 10]*1.023e6);
            if rates.isKey(freqName)
                cr = rates(freqName);
            else
                cr = NaN;
            end
        end
        
    end
end

