function [fitresult,gof,yfit] = PseudoVoigtFit(x,y)
%% Lorentzian + Gaussian fit
% Fits a Lorentzian + Gaussian curve to your ESR spectrum 
% plots the results and outputs fitting parameters, goodness of  and
% fitted curvefit

%% Find good starting parameters
[~,Nmin]=min(y);[Max,Nmax]=max(y);
Gamma= x(Nmin)-x(Nmax);
B0=x(Nmax)+Gamma/2;

%% Fit:
[xData, yData] = prepareCurveData( x, y );

% Set up fittype and options.
ft = fittype( '-a*(x-x0)/GammaG^2*exp(-(x-x0)^2/(2*GammaG^2))-b*2*GammaL^2*(x-x0)/(GammaL^2+(x-x0)^2)^2', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Robust = 'LAR';
opts.StartPoint = [Gamma Gamma 2*Max 2*Max B0];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Create a figure for the plots.
figure( 'Name', 'Pseudo-Voigt Fit' );

% Plot fit with data.
subplot( 2, 1, 1 );
h = plot( fitresult, xData, yData );
legend( h, 'Experimental (normalised)', 'Pseudo-Voigt Fit', 'Location', 'NorthEast' );
% Label axes
xlabel x
ylabel y
grid on

% Plot residuals.
subplot( 2, 1, 2 );
h = plot( fitresult, xData, yData, 'residuals' );
legend( h, 'residuals', 'Zero Line', 'Location', 'NorthEast' );
% Label axes
xlabel x
ylabel y
grid on

yfit=feval(fitresult,x);


