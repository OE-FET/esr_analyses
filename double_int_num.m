function [IArea, yInt1, yCorr, yInt2] = double_int_num(x, y, varargin)
%DOUBLE_INT_NUM Numerical double integration of signal.
%
%   Performs a numerical double integration of a given signal. It offers
%   optional base-line correction and smoothing steps prior to numerical
%   integration.
%
%   SYNTAX:
%   [IArea, Int1, yCorr, Int2] = DOUBLE_INT_NUM(x, y)
%   [IArea, Int1, yCorr, Int2] = DOUBLE_INT_NUM(x, y, 'baseline', 'n')
%
%   INPUT(S):
%   x - vector with magnetic field
%   y - vector with ESR signal intensity
%
%   OPTIONAL:
%   'baseline' - Perform a base-line correction? 'y' (default) or 'n'
%
%   OUTPUT(S):
%   IArea - double-integrated signal
%   yInt1 - single-integrated signal
%   yCorr - polynomial baseline-corrected ESR signal intensity
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$

import esr_analyses.*
import esr_analyses.utils.*

% default to baseline correction if not specified
baseline = get_kwarg(varargin, 'baseline', 'y');

% perform smoothing upon request
if strcmp(baseline, 'y')
    % create a new figure
    figure();
    % plot spectrum
    stackplot(x, y);

    str = input('Would you like to smooth the spectrum y/[n]?', 's');
        if isequal(str, 'y')
            ysmooth = smooth(y, 'sgolay', 3);
            y = reshape(ysmooth, size(y));
        end
end

dim = size(y);

% numerical integration with trapez algorithm
yInt1 = cumtrapz(x, y);

yCorr = y;

% perform baseline correction upon request
if strcmp(baseline, 'y')
    % plot result from first integration
    stackplot(x, yInt1);
    str = input('Would you like to perform a polynomial base-line correction y/[n]?', 's');
    if strcmp(str, 'y')
        % calculate baseline
        yInt1 = baseline_corr(x, yInt1);
        % differentiate baseline-corrected first integral
        yCorr = diff(yInt1)/(x(2)-x(1));
        % add lost data point to derivative
        yCorr(end+1,:) = y(end,:);
    end
end

% second integration to calculate total area
yInt2 = cumtrapz(x, yInt1);

% plot result for visual check
if strcmp(baseline, 'y')
    stackplot(x, yInt2);
    pause(0.2);
end

% average over last 3% from second integral for total area
IArea = mean(yInt2(round(97*dim(1)/100):end,:));
IArea = IArea';

end