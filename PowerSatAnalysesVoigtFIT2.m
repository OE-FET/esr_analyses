function [argout] = PowerSatAnalysesVoigtFit2(varargin)
%POWERSATANALYSESVOIGTFIT2 Analyses of CW-ESR power saturation measurements
%
% Calculates power saturation curves of integrated intensity and maximum
% intensity. This program fits the ESR signal with two voigt curve
% derivatives and then performs analytical integration.
%
% Advantage: Long tails of the resonance peak are not negletcted.
% Disadvantage: Resonance curve has to be symmetrical and of a voigt-shape.
%
% All spectra are collected in the structure 'argout.ERSIntensity'.
%
% INPUT(S):
% PowerSatAnalysesVoigtFIT()            - opens gui to select data
% ...('signal_path')                    - path to signal data
% ...('signal_path','bg_path')          - path to signal data, path to
%                                         background data
% ...(x,y,pars)                         - magnetic field, intensity,
%                                         parameters
% ...(...,'init', var0)                 - list with starting points for fit
%
% OUTPUT(S):
% argout - structure containing all output data
%
% DEPENDENCIES:
% BrukerRead.m
% normalise_spectrum.m
% subtract_background.m
% ESRVoigtSimulation.m
% gfactor_determination.m
% mw_mean.m
% get_sample_position.m
% Natural Constants
%
% MODIFICATION(S):
% 18-Jun-18
% Ian Jacobs: Modified from PowerSatAnalysesVoigtFit to fit two signals
% and accept initial fit parameters (var0)

%%                              Read data 
%%=========================================================================

try
    [var0, varargin] = get_varargin(varargin, 'init');
    n = nargin - 2;
catch
    n = nargin;
end

switch n
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
end

if strcmp(pars.YTYP, 'IGD')==1
    Pmw = pars.z_axis * 1E-3; % MW Power in W
else
    error('The specified file is not a 2D data file.');
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

if ~exist('var0','var')
    % get initial parameters for fit
    mid  = round(length(Bmw)/2);
    fit1 = pseudo_voigt_fit(x, y(:, mid), 'deriv', 1);

    Area  = abs(fit1.a);                       % peak area
    FWHM_lorentz  = fit1.FWHM_lorentz;         % FWHM of Lorentzian component
    FWHM_gauss  = fit1.FWHM_gauss;             % FWHM of Gaussian component

    T2   = 1/(gmratio * FWHM_lorentz*1E-4);    % starting point for fit, sec
    T1   = 10*T2;                              % starting point for fit, sec
    B0   = fit1.x0;                            % resonance center in Gauss

    % starting points for fit:
    %       Area   B0 T1 T2 Brms       Area   B0 T1 T2 Brms
    var0 = [Area*0.05 B0 T1 T2 FWHM_gauss Area*0.05 B0 T1 T2 FWHM_gauss];
end

%%                          Perform voigt fit
%%=========================================================================

% grid data for fitting algorithm
[X, Y]  = meshgrid(x, Bmw);
Z       = y; 

% function to minimize: sum of squared errors
fitfunc1 = @(var) abs(var(1))*ESRVoigtSimulation(X , var(2), abs(var(3)), abs(var(4)), Y, abs(var(5)), 1, pars.B0MA*1e4);
fitfunc2 = @(var) abs(var(6))*ESRVoigtSimulation(X , var(7), abs(var(8)), abs(var(9)), Y, abs(var(10)), 1, pars.B0MA*1e4);

fitfunc = @(var) fitfunc1(var) + fitfunc2(var);

sumofsquares = @(var) sum(sum(abs(fitfunc(var) - Z).^2));

% Fit model to data with fminsearch (Nelder Mead algorithm, much better
% convergance than Levenberg Marquard or trust Region)
opt = optimset('TolFun', 1e-12, 'TolX', 1e-12, 'PlotFcns', ...
    @optimplotfval, 'MaxFunEvals', 1e12, 'MaxIter', 1e12);
[ft_rslt, sumofsquares_error] = fminsearch(sumofsquares, var0, opt);

A    = [abs(ft_rslt(1)), abs(ft_rslt(6))];
B0   = [abs(ft_rslt(2)), abs(ft_rslt(7))];
T1   = [abs(ft_rslt(3)), abs(ft_rslt(8))];
T2   = [abs(ft_rslt(4)), abs(ft_rslt(9))];
Brms = [abs(ft_rslt(5)), abs(ft_rslt(10))];

% get best fit curve (in a higher resolution version)
BmwPlot         = linspace(min(Bmw), max(Bmw), 4*length(Bmw));
[XPlot, YPlot]  = meshgrid(x, BmwPlot);
fitfuncPlot     = @(var) abs(var(1))*ESRVoigtSimulation(XPlot, var(2), abs(var(3)), abs(var(4)), YPlot, abs(var(5)), 1, pars.B0MA*1e4) + abs(var(6))*ESRVoigtSimulation(XPlot, var(7), abs(var(8)), abs(var(9)), YPlot, abs(var(10)), 1, pars.B0MA*1e4);
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

f1 = figure('Name', '3D voigt Fit');
figure(f1);
hold off;
scatter3(Xt(:), Yt(:), Z(:), '.k');
hold on;
surf(XPlot, YPlot, zFit','FaceAlpha', 0.5, 'EdgeColor','none');
legend('Data', 'Fit', 'Location', 'northeast')
xlabel('Magnetic Field [G]')
ylabel('Microwave Magnetic Field [T]')
zlabel('ESR signal [a.u.]')

f2 = figure('Name', '3D voigt Fit, x-section');
figure(f2);
offset = max(max(y)) * 0.8;
hold off;
stackplot(x, y, 'yoffset', offset, 'style', '.k');
hold on;
stackplot(x, fitfunc(ft_rslt), 'yoffset', offset, 'style','-r');
stackplot(x, fitfunc1(ft_rslt), 'yoffset', offset, 'style','--g');
stackplot(x, fitfunc2(ft_rslt), 'yoffset', offset, 'style','--b');
legend('Data');

%%                      Susceptibility Calculation
%%=========================================================================
for i=1:2
    FWHMLorentz(:,i) = 2/(gmratio*T2(i)) * sqrt(1 + gmratio^2*Bmw.^2*T1(i)*T2(i)); % in Tesla
    s = 1246.572123192064; % scaling factor for pseudo modulation
    FittedAreas(:,i) = s * pars.B0MA * 1e4 * Bmw./FWHMLorentz(:,i) .* A(i);
    GFactor(i) = b2g(B0(i) * 1e-4, pars.MWFQ);
    pars.GFactor = GFactor(i);
    
    [Chi(:,i), dChi(:,i)]     = susceptebility_calc(FittedAreas(:,i), pars);
    [NSpin(:,i), dNSpin(:,i)] = spincounting(FittedAreas(:,i), pars);
end
pars.GFactor = GFactor;

%%                                Output
%%=========================================================================

% create output structure
argout.T        = str2double(strtrim(regexprep(pars.Temperature, 'K', '')));

argout.x        = x;
argout.y        = y;
argout.yFit     = fitfunc(ft_rslt);
argout.yFit1    = fitfunc1(ft_rslt);
argout.yFit3    = fitfunc2(ft_rslt);
argout.Pars     = pars;

argout.A        = A;
argout.T1       = T1;
argout.T2       = T2;
argout.Brms     = Brms;

argout.dA       = [abs(se(1)), abs(se(6))];
argout.dB0      = [abs(se(2)), abs(se(7))];
argout.dT1      = [abs(se(3)), abs(se(8))];
argout.dT2      = [abs(se(4)), abs(se(9))];
argout.dBrms    = [abs(se(5)), abs(se(10))];


argout.B0       = B0;
argout.gfactor  = pars.GFactor;
argout.Bmw      = Bmw;

argout.Chi      = Chi;
argout.dChi     = dChi;
argout.NSpin    = NSpin;
argout.dNSpin   = dNSpin;

argout.Bp2p     = Bp2p;

end
