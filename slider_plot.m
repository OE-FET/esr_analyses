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

if length(x) ~= length(y)
    error('x and y data must have the same length.')
end

S.x = x; S.y = y;
S.style = get_kwarg(varargin, 'style', '');

S.fh = gcf;
S.ax = gca;

plot(S.x, S.y(:,1), S.style);
xlim([min(S.x) max(S.x)]);
ylim(1.1*[min(S.y, [], 'all') max(S.y, [], 'all')]);
grid on;

S.sl = uicontrol('style', 'slide', ...
                 'unit', 'pix', ...
                 'min', 1, 'max', size(S.y, 2), 'val', 1, ...
                 'sliderstep', [1/size(S.y, 2) 1/size(S.y, 2)], ...
                 'callback', {@sl_call, S});
end

function [] = sl_call(varargin)
% Callback for the slider.
[h, S] = varargin{[1, 3]};  % calling handle and data structure.
cla

slice_number = round(get(h, 'value'));
plot(S.ax, S.x, S.y(:,slice_number), S.style);
xlim([min(S.x) max(S.x)]);
ylim(1.1*[min(S.y, [], 'all') max(S.y, [], 'all')]);
grid on;

end