function [fitresult] = saturation_area_fit(Bmw, DoubleIntI, g)

%% Power saturation fit
%
% Fits power saturation of double integrated intensity to determine T1*T2
% with euqation:
%
% A*x/(gmratio*sqrt(1+gmratio^2*x^2*T1*T2)).
%
% A Markov chain Monte Carlo method is used for fitting.
%
% INPUT(S): 
% Bmw - microwave magnetic field in Tesla
% DoubleIntI - double integrated ESR signal
% g - optional: average (isotropic) g-value of sample. If not given, the
%   free electron g-value is used.
%
% OUPUT(S):
% fitresult - structure containing all fitting parameters
% gof - parameters giving quality of fit
%
%%

x = Bmw; y = DoubleIntI;

% calculate gyromagnetic ratio used for the fit
gm = gmratio;
if nargin==3
    gm = g*bmagn / hbar;
end

% determine starting values
A0 = mean(1.1*DoubleIntI*gm./Bmw);
T0 = mean( (A0^2*Bmw(end).^2 - DoubleIntI(end).^2*gm^2)./(Bmw(end).^2.*DoubleIntI(end).^2*gm^4) );

var0 = [A0, T0];

% function to minimize: sum of squared errors
fitfunc = @(var) abs(var(1)) * x ./( gm * sqrt(1+gm^2*x.^2*abs(var(2))) );
sumofsquares = @(var) sum(sum( abs(fitfunc(var) - y).^2  ));

% Fit model to data with fminsearch
opt = optimset('TolFun', 1e-9,'TolX', 1e-9,'PlotFcns', @optimplotfval, 'MaxFunEvals', 1e9, 'MaxIter', 1e9);
[ft_rslt, sumofsquares_error] = fminsearch(sumofsquares, var0, opt);

yFit = fitfunc(ft_rslt);

A = abs(ft_rslt(1));
T1T2 = abs(ft_rslt(2));

%% Estimate fitting errors
dof = length(x) - length(var0); % degrees of freedom in fitting problem
sdr = sqrt(sumofsquares_error/dof); % standard deviation of residuals
J = jacobianest(fitfunc, ft_rslt); % jacobian matrix
Sigma = sdr^2*inv(J'*J); % covariance matrix
se = sqrt(diag(Sigma))'; % parameter standrad errors

dA = abs(se(1));
dT1T2 = abs(se(2));

% Plot results
plot(x,y,'.');
hold on
plot(x, yFit);

% create output structure
fitresult.A = A;
fitresult.dA = dA;
fitresult.T1T2 = T1T2;
fitresult.dT1T2 = dT1T2;
end