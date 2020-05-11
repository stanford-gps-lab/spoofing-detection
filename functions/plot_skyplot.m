%% function plot_skyplot(SV_ids)
% 
% Plots a skyplot showing all the satellites identified by SV_ids.
% 
% Uses the global variables reciever_data, settings, plt from
% readublox_nullsteering.m
% 
% Input:
% SV_ids    Indexes of all the satellites that are to be plotted. These
%           indexes identify the columns of the receiver_data variables to
%           be used.


function plot_skyplot(azTrue, elTrue, azMeas, elMeas)

if nargin < 3
    azMeas = NaN(size(azTrue));
end
if nargin < 4
    elMeas = elTrue;
end

fs = 16;

% convert to polar coords
rT = 90 - elTrue * 180/pi;
rM = 90 - elMeas * 180/pi;
satxT = rT.*cos(azTrue); % true x pos
satxM = rM.*cos(azMeas); % measured x pos
satyT = rT.*sin(azTrue); % true y pos
satyM = rM.*sin(azMeas); % measured y pos


polarhg([30, 60]); hold on;
pm = plot(reshape(satxM, numel(satxM), 1), ...
          reshape(satyM, numel(satyM), 1), ...
          'b+', 'MarkerSize', 8, 'LineWidth', 2);
pt = plot(satxT, satyT, 'r.', 'MarkerSize', 25);

legend([pt; pm], 'Ephemeris', 'Measured', ...
    'FontSize', fs, 'Location', 'best', 'Interpreter', 'latex')

% % save image as matlab figure and .png
% if plt_save
%     save_figure('', plt_path, 'skyplot_estimate')
% end

end








