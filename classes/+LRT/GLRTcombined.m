classdef GLRTcombined < LRT.GLRT
    %class of Generalized Likelihood Ratio Test for comination of variables
    %   Calculates a likelihood ratio test of a combination of variables.
    %   The variables can be of simple vs. simple and simple vs. composite
    %   hypothesis but are all expected to be Gaussian.
    %   
    
    
    properties
        SStests         % vector of Gaussian simple vs. simple LRTs
        SCtests         % vector of simple vs. composite GLRTs
        SC_GLRT         % single GLRT combining all SC GLRTs
        SSmu            % joint mean of all Gaussian simple vs. simple LRTs
        SSSigma         % joint Var of all Gaussian simple vs. simple LRTs
%         SCdof           % joint deg. of freedom of simple vs. composite GLRTs
%         N               % totala original number of measurements
%         p_yH0           % conditional probability of nominal hypothesis
%         p_yH1           % conditional probability of spoofed hypothesis
%         logLambda       % joint log Lambda of all tests
        mean            % mean of distribution under H0
        var             % variance of distribution under H0
        std             % standard deviation of distribution under H0
    end
    
    properties (Dependent)
    end
    
    
    methods
        function obj = GLRTcombined(SStests, SCtests)
            %GLRTcombined(SStests, SCtests)
            %   Constructs an instance of the class for a vector of simple
            %   vs. simple LRTs and simple vs. composite LRTs.
            %   SStest need to be LRT.SS, SCtests need to be of 
            %   SC_GLRT class or each subclasses thereof.
            
            % check inputs
            if nargin < 2
                SCtests = [];
                if nargin < 1
                    SStests = [];
                end
            end
            
            % gather SigmaInv variables
            SSSigmaInv = arrayfun(@(x) x.SigmaInv, SStests(:), ...
                'UniformOutput', false);
            SSy = arrayfun(@(x) x.y, SStests(:), ...
                'UniformOutput', false);
            
            SCSigmaInv = arrayfun(@(x) x.SigmaInv, SCtests(:), ...
                'UniformOutput', false);
            SCy = arrayfun(@(x) x.y, SCtests(:), ...
                'UniformOutput', false);
            
            emptyYs = cellfun(@isempty, [SSy; SCy]);
            if any(emptyYs) && ~all(emptyYs)
                % some but not all GLRTs had measurements attached. Add NaN
                % measurements to other GLRTS for get methods to work
                for c = find(emptyYs)
                    if c <= numel(SStests)
                        SSy{c} = NaN(SStests(c).N, 1);
                    else
                        c1 = c-numel(SStests);
                        SCy{c1} = NaN(SCtests(c1).N, 1);
                    end
                end
            end
            
            % call superclass constructor
            obj@LRT.GLRT(blkdiag(SSSigmaInv{:}, SCSigmaInv{:}), ...
                         vertcat(SSy{:}, SCy{:}));
            
            % collect GLRT objects
            obj.SStests  = SStests(:);
            obj.SCtests  = SCtests(:);
            
            % compute SC GLRT representing all passed SC tests
            mu0s = arrayfun(@(x) x.mu0, obj.SCtests, ...
                'UniformOutput', false);
            obj.SC_GLRT = LRT.SC_GLRT(blkdiag(SCSigmaInv{:}), ...
                                      vertcat(mu0s{:}), ...
                                      sum(arrayfun(@(x) x.p, obj.SCtests)));
            
            % compute joint mu, Sigma of simple vs. simple tests
            [SS_mus, SS_Sigs] = arrayfun(@(x) x.nominalDistribution, obj.SStests);
            obj.SSmu = sum(SS_mus);
            obj.SSSigma = sum(SS_Sigs);
            
            
            % compute joint conditional probabilities (now done by
            % superclass get method)
%             obj.p_yH0 = ...
%                   prod(arrayfun(@(x) x.p_yH0, obj.SStests), 'omitnan') ...
%                 * prod(arrayfun(@(x) x.p_yH0, obj.SCtests), 'omitnan');
%             obj.p_yH1 = ...
%                   prod(arrayfun(@(x) x.p_yH1, obj.SStests), 'omitnan') ...
%                 * prod(arrayfun(@(x) x.p_yH1, obj.SCtests), 'omitnan');
%             
%             % compute joint log Lambda
%             obj.logLambda = ...
%                   sum(arrayfun(@(x) x.logLambda, obj.SStests), 'omitnan') ...
%                 + sum(arrayfun(@(x) x.logLambda, obj.SCtests), 'omitnan');
            
%             % compute total number of measurements
%             obj.N = sum(arrayfun(@(x) x.N, obj.SStests)) ...
%                 + sum(arrayfun(@(x) x.N, obj.SCtests));
                      
            % compute mean, var, std of the distribution under H0
            obj.mean = sum([obj.SSmu,  - obj.SC_GLRT.dof/2], 'omitnan');
            obj.var  = sum([obj.SSSigma, obj.SC_GLRT.dof/2], 'omitnan');
            obj.std  = sqrt(obj.var);
            
        end
        
        % ------------ GET METHODS ---------------------
        % - none -
        
        % ------------ HELPER METHODS -------------------
        
        function xS = inverseCDF(obj, pdf, a, b, p, ep)
            %xS = obj.inverseCDF(pdf, a, b, p, ep)
            %   Calculates the inverse cdf of the pdf.
            %   Simple line search algorithm between a and b. Finds the
            %   point at which the integral of f starting at a is equal to
            %   p. This is equal to the inverse cdf of a distribution
            %   defined by the pdf. Interval of a and b must be chosen to 
            %   contain the solution.
            %   
            %   @params
            %   pdf     <function handle> univariate pdf
            %   a       <double> lower bound on solution
            %   b       <double> upper bound on solution
            %   p       <double> probability for which to solve
            %   ep      <double> (optional) convergence tolerance, 
            %                    default = 1e-6
            %   
            %   @output
            %   xS      <double> quantile, such that int(pdf, -inf, xS) = p

            
            % define optional input: convergence tolerance
            if nargin < 6
                ep = 1e-6;
            end
            
            if length(p) > 1
                xS = zeros(size(p));
                [p, I] = sort(p);
                for i = 1:length(p)
                    if i == 1
                        xS(I(i)) = obj.inverseCDF(pdf, a, b, p(i), ep);
                    else
                        % leverage previous work
                        xS(I(i)) = obj.inverseCDF(pdf, xS(I(i-1)), b, p(i), ep);
                    end
                end
            else
                % set max iterations for convergence while loop
                n_max = 1e3;

                % ensure a <= b
                if a > b
                    bTemp = b;
                    b = a;
                    a = bTemp;
                end

                % calculate integral at lower bound
%                 ya = integral(pdf, -inf, a);
                n = 0;
                while abs(a-b) > ep
                    % evaluate midpoint
                    c = (a+b)/2;
                    yc = integral(pdf, -inf, c, 'RelTol', 1e-10, 'AbsTol', 1e-12);
                    if yc > p % yc is above target, c new upper bound
                        b = c;
                    else % yc below target, c new lower bound
                        a = c;
%                         ya = yc;
                    end
                    % increase counter
                    n = n+1;
                    if n > n_max
                        warning(['Inverse cdf not converged after ', ...
                            num2str(n_max), ' iterations.'])
                    end
                end

                xS = (a + b) / 2;
            end
            
        end
        
        function H0err = H0error(obj, y)
            %obj.H0error(y) calculate the difference between a measurement 
            %y and its mean under nominal conditions.
            %   
            
            H0err = NaN(size(y));
            
            if ~any(arrayfun(@(x) isa(x, 'LRT.SS_periodic'), obj.SStests))
                % otherwise it gets much more complicated...
                
                k = 0;
                % stack H0errors from SS tests
                for GLRT = obj.SStests'
                    if isa(GLRT, 'LRT.SS_periodic')
                        H0err(k+1:k+GLRT.N-1, :) = GLRT.H0error(y(k+1:k+GLRT.N, :));
                        k = k+GLRT.N-1;
                        H0err = H0err(1:end-1, :); % reduce size of H0err
                    else
                        H0err(k+1:k+GLRT.N, :) = GLRT.H0error(y(k+1:k+GLRT.N, :));
                        k = k+GLRT.N;
                    end
                end
                for GLRT = obj.SCtests'
                    H0err(k+1:k+GLRT.N, :) = GLRT.H0error(y(k+1:k+GLRT.N, :));
                    k = k+GLRT.N;
                end
            end
            
            % vectorized:
%             SSerrors = arrayfun(@(x) x.H0error(3*ones(x.N, 1)), ...
%                 obj.SStests, 'UniformOutput', false);
%             SCerrors = arrayfun(@(x) x.H0error(3*ones(x.N, 1)), ...
%                 obj.SStests, 'UniformOutput', false);
%             
%             H0err = [vertcat(SSerrors{:}); vertcat(SCerrors{:})];
            % problem: how to access correct y elements?
%             mu0s = arrayfun(@(x) x.mu0, [obj.SStests; obj.SCtests], ...
%                 'UniformOutput', false);
%             H0err = y - vertcat(mu0s{:});
            % this does not use the classes H0error function
        end
        
        function H1err = H1error(obj, y)
            %obj.H1error(y) calculate the difference between a measurement 
            %y and its mean under nominal conditions.
            %   
            
            H1err = NaN(obj.N, 1);
            
            if ~any(arrayfun(@(x) isa(x, 'LRT.SS_periodic'), obj.SStests))
                % otherwise it gets much more complicated...
                
                k = 0;
                % stack H1errors from SS tests
                for GLRT = obj.SStests'
                    if isa(GLRT, 'LRT.SS_periodic')
                        % output will be one less than number of measurements
                        H1err(k+1:k+GLRT.N-1) = GLRT.H1error(y(k+1:k+GLRT.N));
                        k = k+GLRT.N-1;
                        H1err = H1err(1:end-1);
                    else
                        H1err(k+1:k+GLRT.N) = GLRT.H1error(y(k+1:k+GLRT.N));
                        k = k+GLRT.N;
                    end
                end
                for GLRT = obj.SCtests'
                    H1err(k+1:k+GLRT.N) = GLRT.H1error(y(k+1:k+GLRT.N));
                    k = k+GLRT.N;
                end
            end
        end
        
        function p = getP_logLambdaH0(obj, z)
            %obj.getP_logLambdaH0(z) calculate p(log Lambda | H_0)
            %   Evaluates the joint pdf of all tests at z. z can be a
            %   scalar or vector.
            sig = sqrt(obj.SSSigma);
            k = obj.SC_GLRT.dof;
            
            if ~isfinite(sig) || isempty(obj.SStests) % no normally distributed variables
                p = obj.SC_GLRT.getP_logLambdaH0(z);
            else
                zh = (z(:) - obj.SSmu) / sig;
                if k == 0 % no chi2 distributed variables
                    p = normpdf(zh) / sig;
                else
%                 l = 0:1:300;
% 
%                 evenLogSummands = (2*l.*log(sqrt(2)*(abs(-zh - sig))) ...
%                     + gammaln(k/4+l) - gammaln(2*l+1) - zh.^2/2);
%                 log0i = -zh == sig;
%                 if any(log0i)
%                     evenLogSummands(log0i, :) = ...
%                         [gammaln(k/4)-zh(log0i).^2/2, ...
%                          NaN(1, size(evenLogSummands, 2)-1)];
%                 end
%                 oneAexp = 1  + sign(-zh-sig) ...
%                     .* exp(log(2*pi)/2 + log(abs(-zh - sig) ) - log(2*l+1) ...
%                     - betaln(k/4+l, 1/2));
%                 p = sig.^(k/2-1) ./ gamma(k/2) .* 2.^(k/4-1) ./ sqrt(2*pi) ...
%                     .* nansum(exp(evenLogSummands) .* oneAexp, 2);
%                 
                    % alternative: numerical integration (works better in
                    % matlab)
                    cf = 1 / sqrt(2*pi) / sig / 2.^(k/2) / gamma(k/2);
                    p = cf * integral(@(y) ...
                        y.^(k/2-1) .* exp( -((zh+y/2/sig).^2 + y) /2 ), ...
                        0, inf, 'ArrayValued', true);
                end
            end
        end
        
        function p = getP_logLambdaH1(obj, z, varargin)
            %obj.getP_logLambdaH1(z, lambda) calculate p(log Lambda | H_1)
            %   Evaluates the joint pdf of all tests at z for a given
            %   expected noncentrality parameter lambda.
            %   z can be a scalar or vector.
            
            sig = sqrt(obj.SSSigma);
            
            if obj.SC_GLRT.dof == 0 % no chi2 distributed variables
                p = normpdf(z(:), obj.SSmu, sig);
            elseif ~isfinite(sig) || isempty(obj.SStests) % no normally distributed variables
                p = obj.SC_GLRT.getP_logLambdaH1(z, varargin{:});
            else
                % numerically integrate convolution integral
                if isempty(varargin)
                    lambda_nc = NaN;
                else
                    lambda_nc = varargin{1};
                end
                fX1 = @(y) normpdf(z+y/2, -obj.SSmu, sig);
                fYnc = @(y) ncx2pdf(y, obj.SC_GLRT.dof, lambda_nc);

                p = integral(@(y) fX1(y) .* fYnc(y), 0, inf, ...
                    'ArrayValued', true);
            end
            
        end
        
        function ga = threshold(obj, P_FAmax, varargin)
            %threshold(obj, P_FAmax, dof) calculates the alarm threshold
            %   gamma = obj.threshold(P_FAmax, dof)
            %   Calculates the Executive Monitor (EM) alarm threshold that
            %   satisfies the specified false alert probability.
            %   Calculates the threshold through numerical integration of
            %   the joint pdf of log Lambda.
            %   If the object is made up of only SC GLRT objects, this
            %   method can be called with the additional dof argument.
            
            % sort P_FAmax values in ascending order
%             [P_FAmax, I] = sort(P_FAmax);
%             gamma = zeros(size(P_FAmax));
            
            if obj.SC_GLRT.dof == 0 % no chi2 distributed var
                ga = norminv(P_FAmax, obj.SSmu, sqrt(obj.SSSigma));
            elseif ~isfinite(obj.SSSigma) || isempty(obj.SStests) % no normally distributed variables
                ga = obj.SC_GLRT.threshold(P_FAmax, varargin{:});
            else
                % set integration limits
                qNormal = norminv(P_FAmax, obj.SSmu, sqrt(obj.SSSigma));
                lowerLim = -1/2*chi2inv(1-P_FAmax, obj.SC_GLRT.dof) + qNormal;
                
                % line search within [intLim, normalThresh] for threshold
                ga = obj.inverseCDF(@(x) obj.getP_logLambdaH0(x)', ...
                    min(lowerLim), max(qNormal), P_FAmax);
                
            end
        end
        
        function beta = power(obj, P_FAmax, lambda)
            %obj.power(P_FAmax, lambda) calculates the GLRT test power
            %   p = obj.power(P_FAmax, lambda)
            %   Calcualtes the power of the GLRT to detect an attack of
            %   signals for certain expected sum of squares under H1.
            
            % integrate pdf of distribution under H1 until the threshold
            gamma = obj.threshold(P_FAmax);
            fHm = -obj.SSmu - obj.SC_GLRT.dof/2 - lambda/2; % mean of f_z | H_1
            fHs = (obj.SC_GLRT.dof+2*lambda)/2 + obj.SSSigma; % var of f_z | H_1
            
            if obj.SC_GLRT.dof == 0 % no chi2 distributed var
                beta = normcdf(sqrt(obj.SSSigma) + norminv(P_FAmax));
            elseif ~isfinite(obj.SSSigma) || isempty(obj.SStests) % no normally distributed variables
                beta = obj.SC_GLRT.power(P_FAmax, lambda);
            
            elseif fHm + 10*sqrt(fHs) < gamma
                beta = 1; % approximately right
            else
%                 % get pdf under H1
%                 fH1 = @(z) obj.getP_logLambdaH1(z, lambda);
%                 % calculate cdf value at the treshold
%                 beta = integral(fH1, -inf, gamma, 'ArrayValued', true);
                
                % solve double integral (faster than two individual
                % integrals)
                fNorm = @(x, y) normpdf(y+x/2, -obj.SSmu, sqrt(obj.SSSigma));
                fChi2 = @(x) ncx2pdf(x, obj.SC_GLRT.dof, lambda);
                beta = integral2(@(x, y) fNorm(x, y) .* fChi2(x), 0, inf, -inf, gamma);
            end
        end
        
        % -------- PLOTTING METHODS -------------------
        
        function f = plotLambdaDistribution(obj, varargin)
            %f = obj.plotLambdaDistribution(P_FAmax)
            %   Plot the Lambda decision space with nominal and spoofed
            %   distribtions, thresholds and measured slot.
            
            % use superclass method
            f = plotLambdaDistribution@LRT.GLRT(obj, varargin{:});
            HoldFlag = ishold;
            
            fs = 20;
            
            % add normal and chi2 distribution to the plot
            hold on;
            p0 = f.Children(end);
            if ~isempty(obj.SStests) && ~isempty(obj.SCtests)

                SSsigma = sqrt(obj.SSSigma);
                pdfXvals = min(p0.XData) : 0.01 : obj.SSmu+4*SSsigma;
                p1 = plot(pdfXvals, normpdf(pdfXvals, obj.SSmu, SSsigma), ...
                          'LineWidth', 1.5);
                xp2 = (0 : -0.01 : min(p0.XData));
                p2 = plot(xp2, obj.SC_GLRT.getP_logLambdaH0(xp2), ...
                          'LineWidth', 1.5);
                      
                legend([p0; p1; p2], ...
                    {'Joint pdf'; 'simple vs. simple'; 'simple vs. composite'}, ...
                    'Location', 'best', 'FontSize', fs-3, 'Interpreter', 'latex')
            end
            
            % create title
            titlestr = {'Combined GLRT from '; ...
                [num2str(length(obj.SStests)), ' simple vs. simple and ']; ...
                [num2str(length(obj.SCtests)), ' simple vs. composite tests']};
            
            title(titlestr, ...
                'FontSize', fs-2, 'Interpreter', 'latex')
            
            % reset hold flag
            if ~HoldFlag
                hold off;
            end
        end
        
    end
    
end

