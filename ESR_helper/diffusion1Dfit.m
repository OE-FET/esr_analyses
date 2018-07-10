function [fitresult,gof,yfit] = diffusion1Dfit(x,y)
%% Lorentzian + Gaussian fit
% Fits a Lorentzian + Gaussian curve to your ESR spectrum 
% plots the results and outputs fitting parameters, goodness of  and
% fitted curvefit

%% Find good starting parameters
[~,Nmin]=min(y);[Max,Nmax]=max(y);
Gamma = x(Nmin)-x(Nmax);
B0=x(Nmax)+Gamma/2;

%% Fit:
[xData, yData] = prepareCurveData( x, y );

% Set up fittype and options.
ft = fittype( 'a*diffusion1D(x,x0,width,1)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Robust = 'LAR';
opts.Lower = [0 -Inf -Inf];
opts.StartPoint = [Gamma Gamma*Max B0];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Create a figure for the plots.
figure( 'Name', 'Voigt Fit' );
yfit=feval(fitresult,x);

% Plot fit with data.
figure(1);
subplot( 2, 1, 1 );
plot(x,y,'.',x,yfit,'-');
legend('Experimental (normalised)', '1D difusion model', 'Location', 'SouthWest' );
% Label axes
xlabel 'B (G)'
ylabel 'ESR Intesity'

% Plot residuals.
subplot( 2, 1, 2 );
h = plot(x,y-yfit,'.',x,zeros(1,length(x)),'-');
legend( h, 'residuals', 'Zero Line', 'Location', 'NorthEast' );
% Label axes
xlabel 'B (G)'
ylabel 'Residuals'

display(fitresult);

