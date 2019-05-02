function [output, outputarray] = SliceAnalysesVoigtFit(varargin)
%ESRANALYESFIT performs normalization and spin-counting of an ESR signal by
%fitting it to a voigt function.
%   The ESR signal is normalised according to measurement conditions. If
%   required, a background signal can be subtracted before performing the
%   analyses.
%
%   The total ESR intensity and spin susceptibility are determined by fitting
%   the ESR resonance to a Pseudo-voigt function and subsequent analytical
%   double integration. The number of spins in the sample is then calculated 
%   by assuming that the sample follows a Curie-Weiss law. This assumption is
%   valid when T >> E/2. X-Band EPR typically operates around B = 350 mT while
%   a field of up to 0.5 T corresponds to E/2 = 28.9 ueV. This still remains
%   far lower than the corresponding thermal energy of 450 ueV at 5 K. 
%
%   INPUT(S):
%   ESRAnalysesFIT()            - prompts user for spectrum file
%   ...FIT('Path')              - path to file with ESR data
%   ...FIT('PathSIG','PathBG')  - path to signal, path to background
%   ...FIT(x,y,Pars)            - field, signal, and spectral params
%
%   OUTPUT(S):
%   output                      - output structure containing the
%                                 normalized spectra, measurement conditions,
%                                 fitting parameters with errors, and the
%                                 calculated number of spins and
%                                 susceptibility
%   outputarray                 - array containing the fitting parameters
%                                 and errors for easy copying
%
%   DEPENDENCIES:
%   subtract_background.m
%   normalise_spectrum.m
%   gmarker_calib.m
%   gfactor_determination.m
%   spincounting.m
%   num2clip.m
%   pseudo_voigt_fit.m
%   get_sample_position.m
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2018/07/05 12:58 $    $Revision: 1.1 $

%% Load and normalise spectrum. Subtract a background signal if requested.
close all

switch nargin
    case 0
        str = input('Would you like to subtract a background signal? y/[n]: ','s');
        if strcmpi(str,'y') == 1
            [x, y, Pars] = subtract_background;
        else
            [x, y, Pars] = BrukerRead;
            [x, y, Pars] = normalise_spectrum(x, y, Pars);
        end
    case 1
        Path = varargin{1};
        str = input('Would you like to subtract a background signal? y/[n]: ','s');
        if strcmpi(str,'y') == 1
            [x, y, Pars] = subtract_background(Path);
        else
            [x, y, Pars] = BrukerRead(Path);
            [x, y, Pars] = NormalixseSpectrum(x, y, Pars);
        end
    case 2
        PathSIG = varargin{1}; PathBG = varargin{2};
        [x, y, Pars] = subtract_background(PathSIG, PathBG);
    case 3
        x = varargin{1}; y = varargin{2}; Pars = varargin{3};
        [x, y, Pars] = normalise_spectrum(x, y, Pars);
end

%% Begin fitting

% If a marker is used for g-factor calibration, normalise x-axis according
% to marker position. Does nothing if no marker signal is detected.
[x, y, Pars] = gmarker_calib(x, y, Pars);

% Try to load sample temperature from parameter file. Prompt user for entry 
% if not found.
try
    T = str2double(strtrim(regexprep(Pars.Temperature, 'K', '')));
catch
    T = input('Please give the measurement temperature in Kelvin: ');
    Pars.Temperature = [num2str(T), ' K'];
end

Pars = get_sample_position(Pars);

% Convert MWPW to magnetic field strength in Tesla
% Use Qref = 7500 and c = 2.0 for Nagoya files
Qref = 8355;
c = 2.2;
cCryo = 2.2*1.3; % calibrated with gMarker

f_mean = mw_mean(Pars);

Bmw = f_mean * cCryo * sqrt(Pars.MWPW) * sqrt(Pars.QValue/Qref) * 1E-4; % in Tesla

% Get starting points for fit 
fit1 = pseudo_voigt_fit(x, y, 'deriv', 1);
Area = abs(fit1.a); % area under curve
FWHM_lorentz = fit1.FWHM_lorentz; % Lorentzian linewidth in Gauss
FWHM_gauss = fit1.FWHM_gauss; % Gaussian linewidth in Gauss
T1 = 1e-20; % starting point for fit, sec
T2 = 1/(gmratio * FWHM_lorentz*1E-4); % starting point for fit, sec
B0 = fit1.x0; % resonance center in Gauss
var0 = [Area B0 T2 FWHM_gauss];

% function to minimize: sum of squared errors
fitfunc = @(var) abs(var(1))*ESRVoigtSimulation(x, var(2), T1, abs(var(3)), Bmw, abs(var(4)), 1, Pars.B0MA*1e4)';
sumofsquares = @(var) sum(sum( abs(fitfunc(var) - y).^2  ));

% Fit model to data with fminsearch (Nelder Mead algorithm, much better
% convergance than Levenberg Marquard or trust Region)
options = optimset('TolFun', 1e-8,'TolX', 1e-8,'PlotFcns', @optimplotfval, 'MaxFunEvals', 1e8, 'MaxIter', 1e8);
[opt_val, sumofsquares_error] = fminsearch(sumofsquares, var0, options);

yFit = fitfunc(opt_val);

A = abs(opt_val(1));
B0 = abs(opt_val(2));
T2 = abs(opt_val(3));
Brms = abs(opt_val(4));

%% Estimate fitting errors
dof = length(x) - length(var0); % degrees of freedom in fitting problem
sdr = sqrt(sumofsquares_error/dof); % standard deviation of residuals
J = jacobianest(fitfunc, opt_val); % jacobian matrix
Sigma = sdr^2*inv(J'*J); % covariance matrix
se = sqrt(diag(Sigma))'; % parameter standrad errors

dA = abs(se(1));
dB0 = abs(se(2));
dT2 = abs(se(3));
dBrms = abs(se(4));

%% Get parameters
FWHMLorentz = 2/(gmratio*T2) * sqrt(1 + gmratio^2*Bmw.^2*T1*T2); % in Tesla
dFWHMLorentz = dT2 * 2/(gmratio*T2^2) * sqrt(1 + gmratio^2*Bmw.^2*T1*T2); % in Tesla

s = 1246.572123192064; % scaling factor for pseudo modulation
Area = s * Pars.B0MA * 1e4 * Bmw./FWHMLorentz .* A;
dArea = s * Pars.B0MA * 1e4 * Bmw./FWHMLorentz .* dA;

% count the number of spins, accounting for MW field distribution in cavity
Pars = get_sample_position(Pars);
Pars.GFactor = b2g(B0*1e-4, Pars.MWFQ);

NSpin = spincounting(Area, Pars);
dNSpin = dArea * NSpin/Area;

% calculate suscepebility
Chi = susceptebility_calc(Area, Pars);

% and error margin
dChi = dArea * Chi/Area ;

% compute signal to noise
SNR = snr(yFit,yFit-y);

Bp2p = FWHMLorentz*1e4/sqrt(3);
dBp2p = dFWHMLorentz*1e4/sqrt(3);

% save data to output array
outputnames = {'T', 'gfactor', 'Brms', 'dBrms', 'Bp2p', 'dBp2p', 'Chi','dChi' ,'NSpin','dNSpin','T2', 'dT2'};
outputarray = {T, Pars.GFactor, Brms, dBrms, Bp2p, dBp2p, Chi, dChi, NSpin, dNSpin, T2, dT2};

for k=1:length(outputarray)
    output.(outputnames{k}) = outputarray(k);
end

output.x = x;
output.y = y;
output.Pars = Pars;

output.SNR = SNR;

%Plot data and fit curve
figure; hold on
plot(x, y, '.', 'DisplayName', 'Experiment');
plot(x, yFit, 'DisplayName', 'Fit');
xlabel('Magnetic field [G]')
ylabel('ESR signal [a.u.]')
legend('show')
xlim([x(1), x(end)]);
grid on;
hold off

end