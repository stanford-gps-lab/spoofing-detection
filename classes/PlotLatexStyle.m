classdef PlotLatexStyle < matlab.mixin.Copyable
    %PlotLatexStyle Superclass to add shortcuts for nice latex-like plots.
    
    properties
        fs = 18; % fontsize
        axisLabelArgs = {'FontSize', 18, 'Interpreter', 'latex'};
    end
    
    methods
        
        function l = latexLegend(obj, varargin)
            %l = obj.latexLegend(labels, ___)
            %   Creates a plot legend. Is called the same way as the normal
            %   legend command. Automatically sets the location to 'best',
            %   the FontSize to 16 and the Interpreter to 'latex'.
            l = legend(varargin{:}, 'Location', 'best', ...
                'FontSize', obj.fs-2, 'Interpreter', 'latex');
        end
    end
    
end

