function [IArea, Int1, yCorr, Int2] = double_int_num(x, y, varargin)
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
%   OUTPUT(S):
%   IArea - double integrated signal
%   Int1  - single integrated signal
%   yCorr - polynomial base line corrected ESR signalintensity
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2019/05/06 12:58 $    $Revision: 1.1 $

% default to baseline-correction if not specified
baseline = get_kwarg(varargin, 'baseline', 'y');

% create a new figure
figure();

% perform smoothing upon request
if isequal(baseline, 'y')
    % plot spectrum itself
    stackplot(x, y);

    str = input('Would you like to smooth the spectrum y/[n]?', 's');
        if isequal(str, 'y')
            ysmooth = smooth(y, 'sgolay', 3);
            y = reshape(ysmooth, size(y));
        end
end

dim = size(y);
% preallocate memory for integrated curve
Int1 = zeros(dim);

% numerical integration with trapez-algorithm
for k = 2:dim(1)
    Int1(k,:) = trapz(x(1:k), y(1:k,:));
end

yCorr = y;

% perform baseline correction upon request
if isequal(baseline, 'y')
    % plot result from first integration
    stackplot(x, Int1);
    str = input('Would you like to perform a polynomial base line correction y/[n]?', 's');
    if strcmp(str, 'y') == 1
        % calculate baseline
        Int1 = baseline_corr(x, Int1);
        % differentiate baseline corrected first integral
        yCorr = diff(Int1)/(x(2)-x(1));
        % add lost data point to derivative
        yCorr(end+1,:) = y(end,:);
    end
end

% second integration to calculate totoal Area
Int2 = zeros(dim);
for k = 2:dim(1)
    Int2(k,:) = trapz(x(1:k), Int1(1:k,:));
end

% plot result for visual check
if strcmp(baseline, 'y') == 1
    stackplot(x, Int2);
    pause(0.2);
end

% average over last 3% from second integral for total area
IArea = mean(Int2(round(97*dim(1)/100):end,:));
IArea = IArea';

end