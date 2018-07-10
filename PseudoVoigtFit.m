function [fitresult,gof,yfit] = PseudoVoigtFit(x,y)
%% Lorentzian + Gaussian fit
% Fits a Lorentzian + Gaussian curve to your ESR spectrum 
% plots the results and outputs fitting parameters, goodness of fit and
% fitted cruves.
%
% If y is a matrix of spectral data with the same x-axis, all spetra are
% fitted independently and fitresult is an array with fitresults.

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2018/07/05 12:58 $    $Revision: 0.1 $

%% Find good starting parameters
yfit = zeros(size(y));

for i = 1:size(y,2)

    [~,Nmin]=min(y(:, i));[Max,Nmax]=max(y(:, i));
    Gamma = x(Nmin)-x(Nmax);
    H0 = x(Nmax)+Gamma/2;

    %% Fit:
    [xData, yData] = prepareCurveData( x, y(:, i) );

    % Set up fittype and options.
    ft = fittype( 'PseudoVoigtDeriv(x,a,w,x0,s)' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Robust = 'LAR';
    opts.Lower = [0 0 0 -Inf];
    opts.Upper = [Inf 0 Inf Inf];
    opts.StartPoint = [2*Max 0.5 Gamma H0]; % a, s, w, x0

    % Fit model to data.
    [fitresult{i}, gof{i}] = fit( xData, yData, ft, opts );

    % Create a figure for the plots.
    f1 = figure( 'Name', 'Voigt Fit' );
    yfit(:, i) = feval(fitresult{i},x);

    % Plot fit with data.
    figure(f1);
    subplot( 2, 1, 1 );
    plot(x, y(:, i), '-', x, yfit(:, i), '-');
    legend('Experimental (normalised)', 'Voigt Fit', 'Location', 'SouthWest' );
    % Label axes
    xlabel 'B (Gauss)'
    ylabel 'ESR Intesity'
    axis([min(x) max(x) min(yfit(:, i))*1.1 max(yfit(:, i))*1.1])

    % Plot residuals.
    subplot( 2, 1, 2 );
    h = plot(x, y(:,i)-yfit(:, i), '-', x, zeros(1,length(x)),'-');
    legend( h, 'residuals', 'Zero Line', 'Location', 'NorthEast' );
    % Label axes
    xlabel 'B (Gauss)'
    ylabel 'Residuals'
    axis([min(x) max(x) min(y(:,i)-yfit(:, i))*1.5 max(y(:,i)-yfit(:, i))*1.5])

    display(fitresult{i});
end

if size(fitresult) == 1
    fitresult = fitresult{1};
end

