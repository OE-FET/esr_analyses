function [ycorr, yfit] = baseline_corr(y)
%BASELINE_CORR Performs a baseline correction on the input data.
%
%   [ycorr, yfit] = BASELINE_CORR(y) performs a baseline fit on the input
%   data. The baseline region can be be selected through a GUI, and the
%   baseline is fitted as a spline through the smoothed data. Adjust the
%   number points to smooth over according to the noise in the data.
%
%   OUTPUT(S):
%   ycorr - baseline corrected spectrum
%   yfit - fitted baseline
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2019/05/06 12:58 $    $Revision: 1.2 $

%% process data
% get spectrum size
dim = size(y);
lst = dim(1);
newdimy = [dim(1), prod(dim(2:end))];
y = reshape(y, newdimy);
x = 1:lst;

%% select methods for interpolation
% use spline interpolation for polynomial baseline fit
method  = 'spline';
% average over 10 points for smoothing before fit
avgpts  = round(lst/100); % 1/100 of length of data

%% select area for baseline fit

ok = false;
% re-promt for user input until user confirms good baseline correction
while ~ok
    % plot spectrum
    stackplot(x, real(y));
    fhandle = gcf;

    % set title of plot
    title('Baseline Fit - Select points that belong to baseline');
    % promt user for input of baseline areas
    fprintf(['\n Select the area of the spectrum,', ...
        '\n by indicating points with the curser.', ...
        '\n Press Enter key when done.\n'])
    % open GUI for input, accept only 2 points
    [a, b] = ginput;                                 %#ok
    % round to integer values
    bounds = round(a.');
    bounds = sort(bounds);

    % if points outside x-axis are selected, replace them with axis limits
    bounds(bounds<1) = 1;
    bounds(bounds>lst) = lst;
    bounds = unique(bounds); % delete duplicates

    pts = unique(bounds);

    % smooth curve by averaging over avgpts neighboring points
    npts = length(x);
    pss = zeros(npts, 2);
    pss(:,1) = x - floor(avgpts/2);
    pss(:,2) = pss(:,1) + avgpts;
    pss(pss < 1) = 1;
    pss(pss > lst) = lst;
    yavg = zeros([npts, newdimy(2)]);
    for n = 1:npts
        yavg(n,:) = mean(y(pss(n,1):pss(n,2),:), 1);
    end

    %% perform baseline fit
    yfit = interp1(pts, yavg(pts, :), x, method);
    if size(yfit, 1) == 1
        yfit = shiftdim(yfit, 1);    % make yfit a column if it is a row vector
    end

    % plot for visual confirmation
    hold on;
    yL = get(gca, 'YLim');
    phandle = stackplot(x, real(yfit), 'yoffset', 0.5*max(max(y)));
    set(phandle, 'Color', 'blue');
    ylim(gca, yL);
    for b = bounds
        line([b b], yL, 'Color', 'c');
    end
    hold off;

    %% promt user for confirmation of fit
    set(fhandle, 'Name', 'Baseline Fit - Verify baseline')
    answer = input('Would you like to redo the fit and reselect baseline points? y/[n] ', 's');
    if strcmpi(answer, 'y')
        ok = false;
    else
        ok = true;
    end
end
% close figure if not already done by user
if any(findobj('Type', 'figure') == fhandle)
    close(fhandle)
end
% subtact baseline
ycorr = y - yfit;
% reshape spectrum to original form
ycorr = reshape(ycorr, dim);
yfit = reshape(yfit, dim);

end