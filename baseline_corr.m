function [ycorr, ybaseline] = baseline_corr(x, y, method)
%BASELINE_CORR Performs a baseline correction on the input data.
%
%   [ycorr, ybaseline] = BASELINE_CORR(x, y) performs a baseline fit on the 
%   input data. The baseline region can be be selected through a GUI, and
%   the baseline is fitted as a spline through the smoothed data. Adjust
%   thenumber points to smooth over according to the noise in the data.
%
%   OUTPUT(S):
%   x - x-axis values
%   y - y-axis values
%   method - If 'all', all points in the specified interval will be used to
%   determine the baseline. If method == 'interval', only the endpoints and
%   centre point will be used. Default: 'interval'.
%
%   OUTPUT(S):
%   ycorr - baseline-corrected spectrum
%   ybaseline - fitted baseline
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2019/05/06 12:58 $    $Revision: 1.2 $

%% process data
if nargin < 2
    error('Need both x-data and y-data.');
end

if nargin < 3
    method = 'interval';
end

if length(x)~=length(y)
    error("'x' and 'y' must have the same length.");
end

% get spectrum size
N = length(x);

%% select area for baseline fit

% rescale all data for better visibility
scaling = max(abs(y));
y = y./scaling;

% re-promt for user input until user confirms good baseline correction
ok = false;
while ~ok
    % plot spectrum
    fig = figure('Name', 'Baseline fit');
    [~, offsets] = stack_plot(x, real(y)); hold on;
    
    % set title of plot
    title('Baseline Fit - Select intervals that belong to baseline');

    % prompt user for input of baseline areas
    fprintf(['\n Select the area of the spectrum', ...
        '\n by indicating points with the cursor.', ...
        '\n Press the Enter key when done.\n'])
    
    % open GUI for input & accept points until user presses Enter
    mask = false(N, 1);
    xpts = [];
    ypts = [];
    yLim = get(gca, 'YLim');
    while true
        int = zeros(2,1);
        for ii = 1:2
            [xpt, ~, button] = ginput(1);
            if isempty(button)
                break;
            end
            % truncate at x-axis limits
            xpt = max(min(x), xpt); xpt = min(max(x), xpt);
            % plot marker
            line([xpt xpt], yLim, 'Color', 'b');
            int(ii) = xpt;
        end
        if isempty(button)
            break;
        end
        int = sort(int);
        % plot interval
        rectangle('Position', [int(1), yLim(1), diff(int), diff(yLim)], ...
                      'FaceColor', [0 0 1 0.1]);
        % get points in interval
        int_msk = int(1)< x & x < int(2);
        mask = or(mask, int_msk);
        
        % get end and center points
        x_int = x(int_msk); y_int = y(int_msk, :);
        
        xpts(end+1:end+3) = [x_int(1); x_int(round(end/2)); x_int(end)];
        ypts(end+1:end+3, :) = [y_int(1, :); y_int(round(end/2), :); y_int(end, :)];
    end
    
    if ~any(mask)
        ybaseline = y.*scaling;
        ycorr = y.*scaling;
        return
    end
    
    %% perform baseline fit
    if strcmp(method, 'all')
        % average over 10 points for smoothing before fit
        avgpts  = round(N/100); % 1/100 of length of data
        ysmooth = smoothdata(y, 'movmean', avgpts);
        ybaseline = interp1(x(mask), ysmooth(mask, :), x, 'makima');
    elseif strcmp(method, 'interval')
        ybaseline = interp1(xpts, ypts, x, 'spline');
    end

    % make yfit a column if it is a row vector
    if size(ybaseline, 1) == 1
        ybaseline = shiftdim(ybaseline, 1); 
    end

    % plot for visual confirmation
    yL = get(gca, 'YLim');
    hold on;
    phandle = stack_plot(x, ybaseline, 'yoffsets', offsets);
    set(phandle, 'Color', 'blue');
    ylim(gca, yL);

    %% promt user for confirmation of fit
    title('Baseline Fit - Verify the baseline correction');
    verfiyQ = input('Would you like to redo the fit and re-select baseline points? y/[n] ', 's');
    if strcmpi(verfiyQ, 'y')
        ok = false;
    else
        ok = true;
    end
end
% close figure if not already done by user
try                                                                 %#ok
    close(fig)
end

% subtact baseline
ycorr = y - ybaseline;

% undo scaling
ybaseline = ybaseline.*scaling;
ycorr = ycorr.*scaling;

end