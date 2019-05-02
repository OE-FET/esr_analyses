function [argout] = PowerSatAnalysesVoigtFit(varargin)
%POWERSATANALYSESVOIGTFIT Analyses of CW-ESR power saturation measurements
%
%   [argout] = PowerSatAnalysesVoigtFIT(varargin) Calculates power
%   saturation curves of integrated intensity and maximum intensity. This
%   programm fits the ESR signal with a voigt curve derivative and then
%   performs analytical integration.
%
%   Advantage: Long tails of the resonance peak are not negletcted.
%   Disadvantage: Resonance curve has to have a Voigt-shape.
%
%   INPUT(S):
%   PowerSatAnalysesVoigtFIT()            - opens gui to select data
%   ...('signal_path')                    - path to signal data
%   ...('signal_path', 'bg_path')          - path to signal data, path to
%                                           background data
%   ...(x, y, pars)                         - magnetic field, intensity,
%                                           parameters
%
%   OUTPUT(S):
%   argout - structure containing all output data
%
%   DEPENDENCIES:
%   BrukerRead.m
%   normalise_spectrum.m
%   subtract_background.m
%   ESRVoigtSimulation.m
%   gfactor_determination.m
%   mw_mean.m
%   get_sample_position.m
%   Natural Constants
%
% 	MODIFICATION(S):
% 	15-Jan-18
% 	Remington Carey: Added to program description. Changed read data cell 
% 	for clarity. Added two-input case. Changed strcmp to strcmpi for 
% 	case-insensitive selection.

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2018/07/05 12:58 $    $Revision: 1.1 $

%%                              Read data 
%%=========================================================================
switch nargin
case 0
    str = input('Would you like to subtract a background signal? y/[n]: ', 's');
    if strcmpi(str, 'y') == 1
        [x, y, pars] = subtract_background();
    else
        [x, y, pars] = BrukerRead();
    end
case 1
    [x, y, pars] = BrukerRead(varargin{1});
case 2
    PathSIG = varargin{1};
    PathBG  = varargin{2};
    [x, y, pars] = subtract_background(PathSIG, PathBG);
case 3
    x    = varargin{1};
    y    = varargin{2};
    pars = varargin{3};
end

%%                         Calculate MW fields
%%=========================================================================

% Get Q value
try
    QValue = pars.QValue;
catch
    QValue = input('Please give cavity QValue:');
    pars.QValue = QValue;
end

if strcmp(pars.YTYP, 'IGD') == 1
    Pmw = pars.z_axis * 1E-3; % MW Power in W
else
    error('The specified file is not a 2D data file.');
end

% get measurement temperature
if isfield(pars, 'Temperature') == 0
    T = input('Please give measurement temperature in K [default = 298 K]:');
    if isempty(T)
        T = 298;
    end
    pars.Temperature = sprintf('%.1f K', T);
end

% ask for height and length of sample if not in Pars
pars    = get_sample_position(pars);

% convert MWPW to magnetic field strength in Tesla
% use Qref = 7500 and c = 2.0 for Nagoya files
Qref    = 8355;
c       = 2.2;      % without cryostat mounted
cCryo   = 2.2*1.3;  % with crystat, calibrated with gMarker

f_mean  = mw_mean(pars);
Bmw     = f_mean * cCryo * sqrt(Pmw) * sqrt(QValue/Qref) * 1E-4; % in Tesla

% normalise for measurement conditions
[x, y, pars] = normalise_spectrum(x, y, pars);

%%                      Get starting points for fit
%%=========================================================================

% get initial parameters for fit
mid  = round(length(Bmw)/2);
fit1 = pseudo_voigt_fit(x, y(:,mid), 'deriv', 1);

Area  = abs(fit1.a);                       % peak area
FWHM_lorentz  = fit1.FWHM_lorentz;         % FWHM of Lorentzian component
FWHM_gauss  = fit1.FWHM_gauss;             % FWHM of Gaussian component

T2   = 1/(gmratio * FWHM_lorentz*1E-4);    % starting point for fit, sec
T1   = 10*T2;                              % starting point for fit, sec
B0   = fit1.x0;                            % resonance center in Gauss

var0 = [Area*0.05 B0 T1 T2 FWHM_gauss];    %starting points

% grid data for fitting algorithm
[X, Y]  = meshgrid(x, Bmw);
Z       = y; 


%%                          Perform voigt fit
%%=========================================================================

% function to minimize: sum of squared errors
fitfunc = @(var) abs(var(1))*ESRVoigtSimulation(X, abs(var(2)), abs(var(3)), abs(var(4)), Y, abs(var(5)), 1, pars.B0MA*1e4);
sumofsquares = @(var) sum(sum( abs(fitfunc(var) - Z).^2  ));

% Fit model to data with fminsearch (Nelder Mead algorithm, much better
% convergance than Levenberg Marquard or trust Region)
opt = optimset('TolFun', 1e-12, 'TolX', 1e-12, 'PlotFcns', @optimplotfval, 'MaxFunEvals', 1e12, 'MaxIter', 1e12);
[ft_rslt, sumofsquares_error] = fminsearch(sumofsquares, var0, opt);

A    = abs(ft_rslt(1));
B0   = abs(ft_rslt(2));
T1   = abs(ft_rslt(3));
T2   = abs(ft_rslt(4));
Brms = abs(ft_rslt(5));

% get best fit curve (in a higher resolution version)
BmwPlot         = linspace(min(Bmw), max(Bmw), 4*length(Bmw));
[XPlot, YPlot]  = meshgrid(x, BmwPlot);
fitfuncPlot     = @(var) abs(var(1))*ESRVoigtSimulation(XPlot, B0, abs(var(2)), abs(var(3)), YPlot, abs(var(4)), 1, pars.B0MA*1e4);
zFit            = fitfuncPlot(ft_rslt);

% get peak-2-peak linewidth of simulated spectrum
zFit1 = zFit(:,1);
Bp2p  = x(zFit1 == min(zFit1)) - x(zFit1 == max(zFit1));
 
%%                       Estimate fitting errors
%%=========================================================================

dof     = length(x) - length(var0);      % degrees of freedom in fitting problem
sdr     = sqrt(sumofsquares_error/dof);  % standard deviation of residuals
J       = jacobianest(fitfunc, ft_rslt); % jacobian matrix
Sigma   = sdr^2*inv(J'*J);               % covariance matrix
se      = sqrt(diag(Sigma))';            % parameter standrad errors

%%                           Plot results
%%=========================================================================

Xt = X'; Yt = Y';

f1 = figure( 'Name', '3D voigt Fit' );
figure(f1);
hold off;
scatter3(Xt(:)./10, Yt(:)*10^6, Z(:), '.k');
hold on;
surf(XPlot./10, YPlot*10^6, zFit', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
legend('Data', 'Fit', 'Location', 'northeast')
xlabel('Magnetic field [mT]')
ylabel('Microwave field [µT]')
zlabel('ESR signal [a.u.]')
axis tight

f2 = figure( 'Name', '3D voigt Fit, x-section' );
figure(f2);
offset = max(max(y))*0.8;
hold off;
stackplot(x, y, 'yoffset', offset, 'style', '.k');
hold on;
stackplot(x, fitfunc(ft_rslt), 'yoffset', offset, 'style', '-r');
legend('Data');

%%                      Susceptibility Calculation
%%=========================================================================
FWHMLorentz = 2/(gmratio*T2) * sqrt(1 + gmratio^2*Bmw.^2*T1*T2); % in Tesla
s = 1246.572123192064;  % scaling factor for pseudo modulation
FittedAreas = s * pars.B0MA * 1e4 * Bmw./FWHMLorentz .* A;
pars.GFactor = b2g(B0*1e-4, pars.MWFQ);

[Chi, dChi]     = susceptebility_calc(FittedAreas, pars);
[NSpin, dNSpin] = spincounting(FittedAreas, pars);

%%                                Output
%%=========================================================================

% create output structure
argout.T        = str2double(strtrim(regexprep(pars.Temperature, 'K', '')));

argout.x        = x;
argout.y        = y;
argout.yFit     = zFit;
argout.Pars     = pars;

argout.A        = A;
argout.T1       = T1;
argout.T2       = T2;
argout.Brms     = Brms;

argout.dA       = se(1);
argout.dT1      = se(2);
argout.dT2      = se(3);
argout.dBrms    = se(4);

argout.B0       = B0;
argout.gfactor  = pars.GFactor;
argout.Bmw      = Bmw;

argout.Chi      = Chi(1);
argout.dChi     = dChi(1);
argout.NSpin    = NSpin(1);
argout.dNSpin   = dNSpin(1);

argout.Bp2p     = Bp2p;

end
