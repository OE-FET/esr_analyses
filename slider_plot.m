function [S] = slider_plot(x, y, varargin)
%SLIDER_PLOT plot different slices of y-data according to slider location.
%
%   INPUT(S):
%   x, y - Data to plot. Arrays must have the same length.
%   'style' - String specifying plot style.
%
%   OUTPUT(S):
%   S - Structure holding references to figure and axis handles as well as
%       the plotted data.
%

import esr_analyses.*
import esr_analyses.utils.*

if size(x,1) ~= size(y,1)
    error('x and y data must have the same length.')
end

S.x = x; S.y = y;
S.style = get_kwarg(varargin, 'style', '');
S.slider_axis = get_kwarg(varargin, 'slider_axis', []);
S.slider_unit = get_kwarg(varargin, 'slider_unit', []);

S.fh = gcf;
S.ax = gca;

plot(S.x, squeeze(S.y(:,1,:)), S.style);
if ~isempty(S.slider_axis)
    dim = [0.2 0 0.3 0.3];
    str = [num2str(S.slider_axis(1)), S.slider_unit];
    S.a = annotation('textbox', dim, 'String', str, 'FitBoxToText','on');
    S.a.EdgeColor = 'none';
end
xlim([min(S.x) max(S.x)]);
ylim(1.1*[min(S.y, [], 'all') max(S.y, [], 'all')]);
grid on;

n_plots = size(S.y, 2);

minor_step_percent = 1/(n_plots - 1);
major_step_percent = 1/(n_plots - 1);

S.sl = uicontrol('style', 'slide', ...
                 'unit', 'pix', ...
                 'min', 1, 'max', n_plots, 'val', 1, ...
                 'sliderstep', [minor_step_percent major_step_percent], ...
                 'callback', {@sl_call, S});
end

function [] = sl_call(varargin)
% Callback for the slider.
[h, S] = varargin{[1, 3]};  % calling handle and data structure.
cla

slice_number = round(get(h, 'value'));

plot(S.ax, S.x, squeeze(S.y(:,slice_number,:)), S.style);
if ~isempty(S.slider_axis)
    str = [num2str(S.slider_axis(slice_number)), S.slider_unit];
    S.a.String = str;
end
xlim([min(S.x) max(S.x)]);
ylim(1.1*[min(S.y, [], 'all') max(S.y, [], 'all')]);
grid on;

end