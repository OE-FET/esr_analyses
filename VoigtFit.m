function [fitresult,gof,yfit] = VoigtFit(x, y, varargin)
%% Lorentzian + Gaussian fit
% Fits a Lorentzian + Gaussian curve to your ESR spectrum 
% plots the results and outputs fitting parameters, goodness of  and
% fitted curvefit

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2018/07/05 12:58 $    $Revision: 0.1 $

%% Find good starting parameters

[~,Nmin] = min(y); [Max,Nmax]=max(y);

try; Gamma = getVarargin(varargin, 'Gamma'); catch; Gamma = x(Nmin)-x(Nmax); end
try; B0 = getVarargin(varargin, 'B0'); catch; B0 = x(Nmax)+Gamma/2; end

%% Fit:
[xData, yData] = prepareCurveData(x, y);

% Set up fittype and options.
ft = fittype('a*Voigt(x,x0,FWHMGauss,FWHMLorentz,1)', 'independent', 'x', 'dependent', 'y');
opts = fitoptions('Method', 'NonlinearLeastSquares');
opts.Display = 'Off';
opts.Robust = 'LAR';
opts.Lower = [0 0 -Inf -Inf];
opts.StartPoint = [Gamma Gamma Gamma*Max B0];

% Fit model to data.
[fitresult, gof] = fit(xData, yData, ft, opts);

% Create a figure for the plots.
f1 = figure('Name', 'Voigt Fit');
yfit = feval(fitresult,x);

% Plot fit with data.
figure(f1);
subplot(2, 1, 1);
plot(x,y,'.',x,yfit,'-');
legend('Experimental (normalised)', 'Voigt Fit', 'Location', 'SouthWest' );
% Label axes
xlabel 'B (Gauss)'
ylabel 'ESR Intesity'

% Plot residuals.
subplot(2, 1, 2);
h = plot(x,y-yfit,'.',x,zeros(1,length(x)),'-');
legend(h, 'residuals', 'Zero Line', 'Location', 'NorthEast');
% Label axes
xlabel 'B (Gauss)'
ylabel 'Residuals'

display(fitresult);

