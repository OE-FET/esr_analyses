function [IArea,Int1,yCorr,Int2] = DoubleIntNUM(x,y,baseline)
%DOUBLEINTNUM Numerical double integration of ESR signal
%   [IArea,Int1,yCorr,Int2] = DOUBLEINTNUM(x,y,baseline) performs a numerical
%   double integration of a given derivative ESR signal. It offers base-line
%   corrections of both the derivative and integrated signal.
%
%   INPUT:
%   x - vector with magnetic field
%   y - vector with ESR signal intensity
%   baseline - Optional. Enter as string 'y' or 'n'. If 'n', the user will not be
%   prompted to perfrom a baseline correction. If no input is given, 'y' is
%   assumed.
%
%   OUTPUT:
%   IArea - double integrated signal
%   Int1 - single integrated signal
%   yCorr - polynomial base line corrected ESR signalintensity
%
%   DEPENDENCIES:
%   StackPlot.m
%   BaselineFit.m
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2018/07/05 12:58 $    $Revision: 1.1 $

% default to baseline-correction if not specified
if nargin < 3
    baseline = 'y';
end

% perform smoothing upon request
if isequal(baseline,'y')   
    % plot spectrum itself
    hold off;
    StackPlot(x,y);

    str = input('Would you like to smooth the spectrum y/[n]?','s');

        if isequal(str,'y')
            ysmooth = smooth(y,'sgolay',3);
            y = reshape(ysmooth,size(y));
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
if isequal(baseline,'y')
    % plot result from first integration
    StackPlot(x,Int1);
    str = input('Would you like to perform a polynomial base line correction y/[n]?','s');
    if strcmp(str,'y') == 1
        % calculate baseline
        Int1 = BaselineFit(Int1);
        % differentiate baseline corrected first integral
        yCorr = diff(Int1)/(x(2)-x(1));
        % add lost data point to derivative
        yCorr(end+1,:) = y(end,:);
    end
end

% second integration to calculate totoal Area
Int2 = zeros(dim);
for k = 2:dim(1)
    Int2(k, :) = trapz(x(1:k), Int1(1:k,:));
end

% plot result for visual check
if strcmp(baseline, 'y') == 1
    StackPlot(x, Int2);
    pause(0.2);
end

% average over last 3% from second integral for total area
IArea = mean(Int2(round(97*dim(1)/100):end,:));
IArea = IArea';

end