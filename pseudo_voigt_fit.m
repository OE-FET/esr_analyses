function [fitresult,gof,yfit] = pseudo_voigt_fit(x, y, varargin)
%PSEUDO_VOIGT_FIT fits a Pseudo-Voigt curve to the given data, plots the 
% results and outputs fitting parameters, goodness of fit and fitted
% cruves.
%
% If y is a matrix of spectral data with the same x-axis, all spetra are
% fitted independently and fitresult is an array with fitresults.

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2018/07/05 12:58 $    $Revision: 0.1 $

deriv = get_varargin(varargin, 'deriv', 0);
yfit = zeros(size(y));

for i = 1:size(y, 2)
    
    [xData, yData] = prepareCurveData(x, y(:, i));
    
    %% Find good starting parameters

    [Min, Nmin] = min(yData); [Max, Nmax] = max(yData);
    
    if deriv == 0
        FWHM_default = (x(end) - x(1))/5;
        x0_default = Nmax;
        Area_default = Max * FWHM_default;
    elseif deriv == 1
        FWHM_default = x(Nmin)-x(Nmax);
        x0_default = x(Nmax)+FWHM_default/2;
        Area_default = 2*(Max - Min)*FWHM_default;   
    end

    FWHM = get_varargin(varargin, 'FWHM', FWHM_default);
    x0 = get_varargin(varargin, 'x0', x0_default);
    Area = get_varargin(varargin, 'Area', Area_default);

    %% Fit:

    % Set up fittype and options.
    ft = fittype(@(a, x0, FWHM_gauss, FWHM_lorentz, x) a*pseudo_voigt(x, x0, FWHM_gauss, FWHM_lorentz, 'deriv', deriv), ...
        'independent', 'x',  'dependent', 'y');
    opts = fitoptions('Method', 'NonlinearLeastSquares');
    opts.Display = 'Off';
    opts.Robust = 'LAR';
    opts.Lower = [0 0 0 0];
    opts.StartPoint = [Area x0  FWHM FWHM];

    % Fit model to data.
    [fitresult{i}, gof{i}] = fit(xData, yData, ft, opts);

    % Create a figure for the plots.
    f1 = figure('Name', 'Voigt Fit');
    yfit = feval(fitresult{i}, xData);

    % Plot fit with data.
    figure(f1);
    subplot(2, 1, 1);
    plot(xData, yData, '.', x, yfit, '-');
    legend('Experimental (normalised)', 'Voigt Fit', 'Location', 'SouthWest' );

    % Plot residuals.
    subplot(2, 1, 2);
    h = plot(xData, yData-yfit, '.', xData, zeros(1, length(x)), '-');
    legend(h, 'Residuals', 'Zero Line', 'Location', 'NorthEast');

    display(fitresult{i});

end

if size(fitresult) == 1
    fitresult = fitresult{1};
    gof = gof{1};
end

end
