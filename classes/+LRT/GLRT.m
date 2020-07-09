classdef (Abstract) GLRT < matlab.mixin.Heterogeneous
    %Abstract class of Generalized Likelihood Ratio Test
    %   Generally used as superclass for any other likelihood ratio test
    %   class. Subclass of matlab.mixin.Hererogenous to allow for
    %   concatination.
    
    properties
        y               % measurement
        SigmaInv        % inverse of measurement covariance matrix
        Phalf           % Choleski decomposition of SigmaInv
        SVnames         % names of involved satellites
        N               % original number of measurements
        consideredSats  % logical indices: sat considered for spoof detection
        excludedSats    % logical indices: sat excluded from consideration
    end
    
    properties (Dependent)
        p_yH0           % conditional probability of nominal hypothesis
        p_yH1           % conditional probability of spoofed hypothesis
        logLambda       % LR test metric = log( p(y|H0) / p(y|H1) )
        sats2use        % logical indices: considered - excluded
    end
    
    methods
        function obj = GLRT(SigmaInv, y, SVnames)
            %GLRT(y) construct a GLRT class
            %   Constructs a GLRT class. Can be called with variable number
            %   of arguments.
            %   
            
            % check if SVnames cell was passed
            if nargin < 3
                SVnames = {};
                % check if measurement vector was passed
                if nargin < 2
                    y = [];
                    % check if inverse covariance matrix passed
                    if nargin < 1
                        SigmaInv = [];
                    end
                end
            end
            
            obj.SigmaInv = SigmaInv;
            obj.y = y;
            obj.SVnames = SVnames;
            
            % start by considering all satellites
            obj.N               = size(SigmaInv, 1);
            obj.excludedSats    = false(obj.N, 1);
            if length(obj.y) == obj.N
                obj.consideredSats = isfinite(obj.y) & obj.y ~= 0;
            else
                obj.consideredSats  = true(obj.N, 1);
            end
            
            [~, S, V] = svd(obj.SigmaInv(obj.consideredSats, obj.consideredSats));
            obj.Phalf = zeros(size(obj.SigmaInv));
            obj.Phalf(obj.consideredSats, obj.consideredSats) = sqrt(S)*V';
            
        end
        
        % ------------ GET METHODS ---------------------
        
        function p_yH0 = get.p_yH0(obj)
            %get.p_yH0 Calculates the conditional probability of y given H0
            %   p_yH0 = get.p_yH0(obj)
            %   Calculates the conditional probability of the measurement
            %   given the nominal Hypothesis.
            %   
            if any(obj.sats2use) && ~isempty(obj.y)
                p_yH0 = obj.getP_yH0(obj.y);
            else
                p_yH0 = NaN;
            end
        end
        
        function p_yH1 = get.p_yH1(obj)
            %get.p_yH1 Calculates the conditional probability of y given H1
            %   p_yH1 = get.p_yH1(obj)
            %   Calculates the conditional probability of the measured
            %   azimuths given the spoofed Hypothesis.
            %   
            if any(obj.sats2use) && ~isempty(obj.y)
                p_yH1 = obj.getP_yH1(obj.y);
            else
                p_yH1 = NaN;
            end
        end
        
        function lL = get.logLambda(obj)
            %get.logLambda Calculates the log likelihoood ratio.
            %   Calculates the Neyman Pearson criteria variable, the log of
            %   the ratio of conditional probabilities p(y|H0) / p(y|H1).
            if any(obj.sats2use) && ~isempty(obj.y)
                lL = obj.getLogLambda(obj.y);
            else
                lL = NaN;
            end
        end
        
        function sats2use = get.sats2use(obj)
            %get.sats2use joins the consideredSats and multipath properties
            %   Calculate the logical indices of the satellites to be used
            %   for probability calculates, considering excluded satellites
            %   and multipath.
            sats2use = obj.consideredSats & ~obj.excludedSats;
        end
        
        % ------------ STANDARD METHODS -------------------
        
        function p_yH0 = getP_yH0(obj, y)
            %getP_yH0 Calculates the conditional probability of y given H0
            %   p_yH0 = obj.p_yH0(y)
            %   Calculates the conditional probability of any measurement
            %   vector given the nominal Hypothesis.
            %   
            s2u = obj.sats2use;
            if any(s2u)
                p_yH0 = sqrt(det(obj.SigmaInv)) ...
                    * mvnpdf((obj.Phalf * obj.H0error(y))');
            else
                p_yH0 = NaN;
            end
        end
        
        function p_yH1 = getP_yH1(obj, y)
            %getP_yH1 Calculates the conditional probability of y given H1
            %   p_yH1 = obj.p_yH0(y)
            %   Calculates the conditional probability of any measurement
            %   vector given the spoofed Hypothesis.
            %   
            s2u = obj.sats2use;
            if any(s2u)
                p_yH1 = sqrt(det(obj.SigmaInv)) ...
                    * mvnpdf((obj.Phalf * obj.H1error(y))');
            else
                p_yH1 = NaN;
            end
        end
        
        function lLy = getLogLambda(obj, y)
            %obj.getLogLambda(y) Computes the log Lambda of a measurement y
            %   For a given measurement y, calculates the resulting value
            %   of the decision variable log Lambda. If no argument is
            %   passed, returns the saved value of obj.logLambda.
            %   y can be a column vector of a matrix of measurements. If a
            %   matrix is passed, considers every column as a measurement.
            
            if nargin < 1
                y = obj.y;
            end
            
            % remove NaN values from sum by setting to zero
            allnans = all(isnan(y), 1); % to later set result to NaN
            y(isnan(y)) = 0;
            
            % vectorized for a matrix of column vector measurements
            lLy = - sum((obj.Phalf * obj.H0error(y)).^2, 1) / 2 ...
                  + sum((obj.Phalf * obj.H1error(y)).^2, 1) / 2;
            
            % set to NaN when only NaNs passed
            lLy(allnans) = NaN;
        end
        
        function alarm = test(obj, y, P_FAmax)
            %alarm = obj.test(y, P_FAmax) Tests to see if measurement y 
            %raises alarm given a maximum false alert probability
            %   Compares log Lambda(y) to the detection threshold resulting
            %   from P_FAmax. Returns a logical output that is true if an
            %   alarm is raised.
            %   
            %   y is a matrix of measurements where each column is
            %   considered a measurement.
            %   P_FAmax is a scalar or vector of multiple P_FAmax values.
            %   Returns a matrix of logicals with one row per P_FAmax value
            %   and one column per column in y.
            
            alarm = obj.getLogLambda(y) < obj.threshold(P_FAmax(:));
            
        end
        
        function pMD = missedDetectionP(obj, P_FAmax, varargin)
            %missedDetectionP(P_FAmax, vararagin) calculates the GLRT P_MD
            %   p = obj.missedDetectionP(P_FAmax, varargin)
            %   Calcualtes the probability of missed detection  of the GLRT
            %   to detect an attack.
            %   Depending on the GLRT, additional parameters need to be
            %   passed next to the false alert probability, the same as to
            %   obj.power(P_FAmax, ...)
            %   For a simple vs. composite GLRt, this is the noncentrality 
            %   parameter under H1, lambda.
            
            pMD = 1 - obj.power(P_FAmax, varargin{:});
            
        end
        
        % ------------ PLOTTING METHODS -------------------
        
        function f = plotLambdaDistribution(obj, P_FAmax)
            %f = obj.plotLambdaDistribution(P_FAmax)
            %   Plot the Lambda decision space with nominal and spoofed
            %   distribtions, thresholds and measured slot.
            
            % get P_FAmax input if given
            if nargin == 1
                P_FAmax = 10 .^ [-9; -8; -7; -6];
            end
            pFApot = log10(P_FAmax);
            gammas = obj.threshold(P_FAmax);

            % check if hold on is active
            HoldFlag = ishold;
            
            
            % set plot parameters
            fs = 20;
            xPlot = min(gammas) : 0.01 : max(20, -min(gammas));
            % plot pdf under H0
            plot(xPlot, obj.getP_logLambdaH0(xPlot), ...
                'LineWidth', 2)
            f = gca;
            hold on; grid on;
            % plot pdf under H1
            plot(xPlot, obj.getP_logLambdaH1(xPlot), ...
                'LineWidth', 2)
            % adjust tick label font size
            f.FontSize = fs-7;
            % draw thresholds
            line([gammas, gammas], f.YLim, ...
                'Color', 'black', 'LineStyle', '--', 'LineWidth', 1.5)

            % plot threshold lines
            for g = 1:length(gammas)
                text(gammas(g), f.YLim(end) - (f.YLim(end)-f.YLim(1))*0.05*g, ...
                    ['$1e', num2str(pFApot(g)), '$'], ...
                    'HorizontalAlignment', 'left', 'Color', 'black', ...
                    'FontSize', fs-4, 'Interpreter', 'latex')
            end

            % plot actual measurement
            line([obj.logLambda, obj.logLambda], f.YLim, ...
                'Color', 'blue', 'LineWidth', 1.5)
            % label plot
            text(obj.logLambda, f.YLim(end) / 2, ...
                    ['$\log\Lambda(y) = ', num2str(obj.logLambda), '$'], ...
                    'HorizontalAlignment', 'right', 'Color', 'blue', ...
                    'FontSize', fs-4, 'Interpreter', 'latex')

            f.XTick = f.XLim(1):10:f.XLim(end);
            f.XTickLabel = cellfun(@(x) num2str(x), num2cell(f.XTick'), ...
                'UniformOutput', false);

            ylabel('pdf of $\log \Lambda(y)$', ...
                'FontSize', fs, 'Interpreter', 'latex')
            xlabel('$\log \Lambda(y)$', ...
                'FontSize', fs, 'Interpreter', 'latex')
%             titlestr = {['y: $', ...
%                 mat2str(round(obj.y(obj.sats2use), 2)), '$']; ...
%                 ['$\sigma_i$ s: $', ...
%                 mat2str(round(1./obj.w(obj.sats2use), 2)), ...
%                 '$']};
            if ~isempty(obj.SVnames)
                title(['Satellites: ', strjoin(obj.SVnames(obj.sats2use), ' ')], ...
                    'FontSize', fs-2, 'Interpreter', 'latex')
            end
            
            % reset hold flag
            if ~HoldFlag
                hold off;
            end
        end
        
        function f = plotMeas(obj, satellites)
            %f = obj.plotMeas(satellites)
            %   Plots the normalized measurement error under nominal
            %   conditions for each satellite.
            %   Defaults to obj.sats2use if no satellite selection is
            %   specified.
            
            if nargin < 2
                satellites = obj.sats2use;
            end
            % check if hold on is active
            HoldFlag = ishold;
            % plot offset from mean
            plot(1:sum(satellites), ...
                obj.Phalf * obj.H0error(obj.y(satellites)), ...
                '+', 'LineWidth', 2)
            f = gca;
            hold on; grid on;
            % plor actual azimuths
            plot(obj.mu0(satellites), 'r.', 'MarkerSize', 25)
            % reset xticks
            f.XTick = 1:sum(satellites);
            if ~isempty(obj.SVnames)
                f.XTickLabel = obj.SVnames(satellites);
            end
            % label
            xlabel('Satellite', 'FontSize', 16, 'Interpreter', 'latex')
            ylabel('Normalized $y - \mu_0$', ...
                'FontSize', 16, 'Interpreter', 'latex')
            legend('Measured, $2\sigma$', 'Actual', ...
                'Location', 'best', 'FontSize', 16, 'Interpreter', 'latex')
            
            if ~HoldFlag
                hold off;
            end
        end
        
    end
    
    % override getDefaultScalarElement as this is an abstract class
    methods (Static, Sealed, Access = protected)
        function default_object = getDefaultScalarElement
            default_object = LRT.SC_GLRT;
        end
    end

    
    % ------------ ABSTRACT METHODS ---------------------
    methods (Abstract)
        
        H0error(obj, y) % calculate the error of y | H0 = y - mu_0
        H1error(obj, y) % calculate the error of y | H1 = y - mu_1
        
        getP_logLambdaH0(obj, logLambda) % evaluate p(logLambda | H0)
        getP_logLambdaH1(obj, logLambda) % evaluate p(logLambda | H1)
        
        threshold(obj, P_FAmax) % calculate the detection threshold
        
        power(obj, P_FAmax) % calculate the detection power
        
    end
    
end

