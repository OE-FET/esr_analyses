function [argout] = PowerSatAnalysesLorentzFit(varargin)
%POWERSATANALYSESLORENTYFIT Analyses of CW-ESR power saturation measurements
%
%   [argout] = POWERSATANALYSESLORENTYFIT(varargin) Calculates power 
%   saturation curves of integrated intensity and maximum
%   intensity. This programm fits the ESR signal with a Lorentzian 
%   derivative and then performs analytical integration.
%
%   Advantage: Long tails of the resonance peak are not negletcted.
%   Disadvantage: Resonance curve has to be a Lorentzian.
%
%   All spectra are collected in the structure 'argout.ERSIntensity'.
%
%   INPUT(S):
%   POWERSATANALYSESLORENTYFIT()          - opens gui to select data
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
%   ESRLorentzSimulation.m
%   gfactor_determination.m
%   mw_mean.m
%   get_sample_position.m
%   Natural Constants
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2018/07/05 12:58 $    $Revision: 0.1 $

%%                              Read data 
%%=========================================================================
switch nargin
case 0
    str = input('Would you like to subtract a background signal? y/[n]: ', 's');
    if strcmpi(str, 'y') == 1
        [x, y, pars] = subtract_background();
    else
        [x, y, pars] = BrukerRead();
        [x, y, pars] = normalise_spectrum(x, y, pars);
    end
case 1
    [x, y, pars] = BrukerRead(varargin{1});
    [x, y, pars] = normalise_spectrum(x, y, pars);
case 2
    PathSIG = varargin{1};
    PathBG  = varargin{2};
    [x, y, pars] = subtract_background(PathSIG, PathBG);
case 3
    x    = varargin{1};
    y    = varargin{2};
    pars = varargin{3};
    [x, y, pars] = normalise_spectrum(x, y, pars);
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

if strcmp(pars.YTYP, 'IGD')==1
    Pmw = pars.z_axis * 1E-3; % MW Power in W
else
    error('The specified file is not a 2D data file.');
end

% get measurement temperature
if isfield(pars, 'Temperature')==0
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
[x, y] = normalise_spectrum(x, y, pars);

%%                      Get starting points for fit
%%=========================================================================

% get initial parameters for fit
fit1 = pseudo_voigt_fit(x, y(:,end));

Ipp = abs(fit1.a); % peak to peak amplitude
Hpp = fit1.w; % peak to peak line width
T2 = 2/sqrt(3) * 1/(gmratio * Hpp*1E-4);
T1 = 2*T2;
B0 = fit1.x0; % resonance center in Gauss
modAmp = pars.B0MA; % p2p modulation amplitude in Tesla

%%                          Perform Lorentz fit
%%=========================================================================

var0 = [Ipp*1e6 T1 T2]; % vector with starting points

[X, Y] = meshgrid(x, Bmw); Z = y; % grid data for fitting algorithm

% function to minimize: sum of squared errors
fitfunc = @(var) var(1)*ESRLorentzSimulation(X, B0, var(2), var(3), Y, 1, modAmp);
sumofsquares = @(var) sum(sum( (fitfunc(var) - Z).^2) );

% Fit model to data with fminsearch (Nelder Mead algorithm, much better
% convergance than Levenberg Marquard or trust Region)
% opt = optimset('TolFun', 1e-7, 'TolX', 1e-7, 'PlotFcns', @optimplotfval);
opt = optimset('TolFun', 1e-7, 'TolX', 1e-7);
[ft_rslt, sumofsquares_error] = fminsearch(sumofsquares, var0, opt);

A = abs(ft_rslt(1));
T1 = abs(ft_rslt(2));
T2 = abs(ft_rslt(3));

% get best fit curve
zFit = ft_rslt(1) * ESRLorentzSimulation(X, B0, ft_rslt(2), ft_rslt(3), Y, 1, modAmp);

% get peak 2 peak linewidth
zFit1 = zFit(:,1);
Bp2p = x(zFit1 == min(zFit1)) - x(zFit1 == max(zFit1));

%% Estimate fitting errors
dof = length(x) - length(var0); % degrees of freedom in fitting problem
sdr = sqrt(sumofsquares_error/dof); % standard deviation of residuals
J = jacobianest(fitfunc, ft_rslt); % jacobian matrix
Sigma = sdr^2*inv(J'*J); % covariance matrix
se = sqrt(diag(Sigma))'; % parameter standrad errors

%% plot results
Xt = X'; Yt = Y';

figure(3);
hold off;
scatter3(Xt(:), Yt(:), Z(:), '.k');
hold on;
surf(X, Y, zFit', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
legend('Data', 'Fit', 'Location', 'northeast')

figure(4);
offset = max(max(y))*0.8;
hold off;
stackplot(x, y, 'yoffset', offset, 'style', '.k');
hold on;
stackplot(x, zFit, 'yoffset', offset, 'style', '-r');
legend('Data');


%% Susceptebility Calculation
s = 1246.572123192064; % scaling factor for pseudo modulation
LorentzArea = A* Bmw * pi ./( T2*gmratio*sqrt(1 + Bmw.^2*T1*T2*gmratio^2) );
FittedAreas = s * pars.B0MA * 1e4 * LorentzArea;
pars.GFactor = b2g(B0*1e-4, pars.MWFQ);

[Chi, dChi] = susceptebility_calc(FittedAreas, pars);
[NSpin, dNSpin] = spincounting(FittedAreas, pars);

%% Output

% determine g-value
g = b2g(B0*1e-4, pars.MWFQ);

% create output structure
argout.x = x;
argout.y = y;
argout.yFit = zFit;
argout.Pars = pars;

argout.A = A;
argout.T1 = T1;
argout.T2 = T2;

argout.dA = se(1);
argout.dT1 = se(2);
argout.dT2 = se(3);

argout.B0 = B0;
argout.gfactor = g;

argout.Chi = Chi(1);
argout.dChi = dChi(1);
argout.NSpin = NSpin(1);
argout.dNSpin = dNSpin(1);

argout.Bp2p = Bp2p;

end

