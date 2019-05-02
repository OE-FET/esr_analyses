function [output, outputarray] = SliceAnalysesVoigtFit2(varargin)
%ESRANALYESFIT2 performs normalization and spin-counting of an ESR signal by
%fitting it to two voigt functions.
%   The ESR signal is normalised according to measurement conditions. If
%   required, a background signal can be subtracted before performing the
%   analyses.
%
%   The total ESR intensity and spin susceptibility are determined by fitting
%   the ESR resonance to two Pseudo-voigt functions and subsequent analytical
%   double integration. The number of spins in the sample is then calculated 
%   by assuming that the sample follows a Curie-Weiss law. This assumption is
%   valid when T >> E/2. X-Band EPR typically operates around B = 350 mT while
%   a field of up to 0.5 T corresponds to E/2 = 28.9 ueV. This still remains
%   far lower than the corresponding thermal energy of 450 ueV at 5 K. 
%
%   INPUT(S):
%   ESRAnalysesFIT2()           - prompts user for spectrum file
%   ...FIT('Path')              - path to file with ESR data
%   ...FIT('PathSIG','PathBG')  - path to signal, path to background
%   ...FIT(x,y,Pars)            - field, signal, and spectral params
%   ...FIT(x,y,Pars,var0)       - field, signal, spectral params, and
%                                 starting fit parameters
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
%
%   Revision: Ian Jacobs <ij255@cam.ac.uk>
%   $Date: 2018/07/18. Added argument for initial fit parameters (var0)

%% load and normalise spectrum, subtract a background spectrum if requested
close all

switch nargin
    case 0
        str = input('Would you like to subtract a background signal? y/[n]','s');
        if strcmp(str,'y')==1
        [x, y, Pars] = subtract_background;
        else
        [x, y, Pars] = normalise_spectrum;
        end
    case 1
        Path = varargin{1};
        str = input('Would you like to subtract a background signal? y/[n]','s');
        if strcmp(str,'y') == 1
            [x, y, Pars] = subtract_background(Path);
        else
            [x, y, Pars] = normalise_spectrum(Path);
        end
        
    case 2
        PathSIG = varargin{1};
        PathBG = varargin{2};
        [x, y, Pars] = subtract_background(PathSIG, PathBG);
    case 3
        x = varargin{1}; y = varargin{2}; Pars = varargin{3};
        [x, y, Pars] = normalise_spectrum(x, y, Pars);
    case 4
        x = varargin{1}; y = varargin{2}; Pars = varargin{3};
        var0 = varargin{4};
        [x, y, Pars] = normalise_spectrum(x, y, Pars);
end

% if a marker is used for g-factor calibration, normalise x-axis according
% to marker position
% this function does nothing if no marker signal is detected
[x, y, Pars] = gmarker_calib(x, y, Pars);

% try to load sample temperature from parameter file
% prompt user for entry if not found
try
    T = str2double(strtrim(regexprep(Pars.Temperature,'K','')));
catch
    T = input('Please give the measurement temperature in Kelvin:');
    Pars.Temperature = [num2str(T), ' K'];
end

Pars = get_sample_position(Pars);

% convert MWPW to magnetic field strength in Tesla
% use Qref = 7500 and c = 2.0 for Nagoya files
Qref = 8355;
c = 2.2;

position_correction = mw_mean(Pars);

Bmw = position_correction * c * sqrt(Pars.MWPW*1e3) * sqrt(Pars.QValue/Qref) * 1E-4; % in Tesla

if ~exist('var0', 'var')
    % Get starting points for fit 
    fit1 = pseudo_voigt_fit(x, y, 'deriv', 1);
    Area = abs(fit1.a); % area under curve
    FWHM_lorentz = fit1.FWHM_lorentz; % peak to peak line width in Gauss
    FWHM_gauss = fit1.FWHM_gauss; % Gaussian linewidth in Gauss
    T1 = 1e-20; % starting point for fit, sec
    T2 = 1/(gmratio * FWHM_lorentz*1E-4); % starting point for fit, sec
    B0 = fit1.x0; % resonance center in Gauss
    var0 = [Area/2 B0 T2  FWHM_gauss Area/2 B0 T2 FWHM_gauss];
end


% function to minimize: sum of squared errors
fitfunc = @(var) abs(var(1))*ESRVoigtSimulation(x, var(2), T1, abs(var(3)), Bmw, abs(var(4)), 1, Pars.B0MA*1e4)' ...
    + abs(var(5))*ESRVoigtSimulation(x, var(6), T1, abs(var(7)), Bmw, abs(var(8)), 1, Pars.B0MA*1e4)';
sumofsquares = @(var) sum(sum( abs(fitfunc(var) - y).^2  ));

% Fit model to data with fminsearch (Nelder Mead algorithm, much better
% convergance than Levenberg Marquard or trust Region)
opt = optimset('TolFun',1e-9,'TolX',1e-9,'PlotFcns',@optimplotfval, 'MaxFunEvals', 1e9, 'MaxIter', 1e9);
[ft_rslt, sumofsquares_error] = fminsearch(sumofsquares, var0, opt);

yFit = fitfunc(ft_rslt);

A = [abs(ft_rslt(1)), abs(ft_rslt(5))];
B0 = [abs(ft_rslt(2)), abs(ft_rslt(6))];
T2 = [abs(ft_rslt(3)), abs(ft_rslt(7))];
Brms = [abs(ft_rslt(4)), abs(ft_rslt(8))];

%% Estimate fitting errorsB
dof = length(x) - length(var0); % degrees of freedom in fitting problem
sdr = sqrt(sumofsquares_error/dof); % standard deviation of residuals
J = jacobianest(fitfunc, ft_rslt); % jacobian matrix
Sigma = sdr^2*inv(J'*J); % covariance matrix
se = sqrt(diag(Sigma))'; % parameter standrad errors

dA = [abs(se(1)), abs(se(5))];
dB0 = [abs(se(2)), abs(se(6))];
dT2 = [abs(se(3)), abs(se(7))];
dBrms = [abs(se(4)), abs(se(8))];

%% Susceptebility Calculation
for i=1:length(var0)/4
    FWHMLorentz(i) = 2/(gmratio*T2(i)) * sqrt(1 + gmratio^2*Bmw.^2*T1*T2(i)); % in Tesla
    dFWHMLorentz(i) = dT2(i).*2/(gmratio*T2(i)^2) * sqrt(1 + gmratio^2*Bmw.^2*T1*T2(i)); % in Tesla
    s = 1246.572123192064; % scaling factor for pseudo modulation
    FittedArea(:,i) = s * Pars.B0MA * 1e4 * Bmw./FWHMLorentz(i) .* A(i);
    dArea(:,i) = s * Pars.B0MA * 1e4 * Bmw./FWHMLorentz(i) .* dA(i);
    Pars.GFactor = b2g(B0(i)*1e-4, Pars.MWFQ);

    [Chi(i), dChi(i)] = susceptebility_calc(FittedArea(:,i), Pars);
    [NSpin(i), dNSpin(i)] = spincounting(FittedArea(:,i), Pars);
end

% and error margin
dChi = dArea .* Chi./FittedArea ;

% compute signal to noise
SNR = snr(yFit,yFit-y);

Bp2p = FWHMLorentz*1e4/sqrt(3);
dBp2p = dFWHMLorentz*1e4/sqrt(3);

% save data to output array
outputarray = {T, Pars.GFactor, Brms, dBrms, Bp2p, dBp2p, sum(Chi), sum(dChi), sum(NSpin), sum(dNSpin), SNR};
outputnames = {'T','gfactor', 'Brms', 'dBrms', 'Bp2p', 'dBp2p', 'Chi','dChi' ,'NSpin','dNSpin','SNR'};

for k=1:length(outputarray)
    output.(outputnames{k}) = outputarray(k);
end

output.x=x;
output.y=y;

output.T2 = T2;

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