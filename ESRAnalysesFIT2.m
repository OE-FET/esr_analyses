function [output, outputarray] = ESRAnalysesFIT2(varargin)
%% performs normalization and spin-counting of ESR signal
% ESR signal is normalised to a MWFQ of 9.6 GHz and according to measuremnt
% conditions.
% The toal ESR intensity and spin suszeptebility are determined by fitting
% the ESR resonance to a Pseudo-Voigt function and subsequent analytical
% double integration. Assuming a Curie-Weiss law, we calculate the number
% of spins in the sample.
% If required, a backgournd signal can be subtracted before performing the
% analyses.
%
% INPUT:
% Path              -   path to file with ESR data
% (x,y,Pars)        -   spectrum data
% (", var0)         -   starting fit parameters
%
% OUTPUT:
% x_norm, y_norm    -   vectors contaning normalised magnetic and ESR signal
%                       data, respectively
% y_base            -   polynomial baseline corrected normalised ESR signal
% NSPin             -   calculated number of spins assuming a Curie-Weiss
%                       suszeptebility
% Chi               -   spin suszeptebility from double integrated intensity
%
% Dependencies:
% SubtractBackground.m
% NormaliseSpectrum.m
% MarkerCalib.m
% gfactor_determination.m
% SpinCounting.m
% num2clip.m
% PseudoVoigtFit.m
% getSamplePosition.m

% MODIFICATION(S):
% 18-Jun-18
% Ian Jacobs: added the ability to specify initial fit parameters (var0)

%% load and normalise spectrum, subtract a background spectrum if requested
close all

switch nargin
    case 0
        str = input('Would you like to subtract a background signal? y/[n]','s');
        if strcmp(str,'y')==1
        [x, y, Pars] = SubtractBackground;
        else
        [x, y, Pars] = NormaliseSpectrum;
        end
    case 1
        Path = varargin{1};
        str = input('Would you like to subtract a background signal? y/[n]','s');
        if strcmp(str,'y') == 1
            [x, y, Pars] = SubtractBackground(Path);
        else
            [x, y, Pars] = NormaliseSpectrum(Path);
        end
        
    case 2
        PathSIG = varargin{1};
        PathBG = varargin{2};
        [x, y, Pars] = SubtractBackground(PathSIG, PathBG);
    case 3
        x = varargin{1}; y = varargin{2}; Pars = varargin{3};
        [x, y, Pars] = NormaliseSpectrum(x, y, Pars);
    case 4
        x = varargin{1}; y = varargin{2}; Pars = varargin{3};
        var0 = varargin{4};
        [x, y, Pars] = NormaliseSpectrum(x, y, Pars);
end

% if a marker is used for g-factor calibration, normalise x-axis according
% to marker position
% this function does nothing if no marker signal is detected
[x, y, Pars] = MarkerCalib(x, y, Pars);

% try to load sample temperature from parameter file
% prompt user for entry if not found
try
    T = str2double(strtrim(regexprep(Pars.Temperature,'K','')));
catch
    T = input('Please give the measurement temperature in Kelvin:');
    Pars.Temperature = [num2str(T), ' K'];
end

Pars = getSamplePosition(Pars);

% convert MWPW to magnetic field strength in Tesla
% use Qref = 7500 and c = 2.0 for Nagoya files
Qref = 8355;
c = 2.2;

position_correction = MWmean(Pars);

Bmw = position_correction * c * sqrt(Pars.MWPW*1e3) * sqrt(Pars.QValue/Qref) * 1E-4; % in Tesla

% get starting points for fit
fit1 = PseudoVoigtFit(x, y);

Ipp = abs(fit1.a); % peak to peak amplitude in Gauss
Hpp = fit1.w; % peak to peak line width in Gauss
T2 = 2/sqrt(3) * 1/(gmratio * Hpp*1E-4); % starting point for fit, sec
T1 = 1e-20; % starting point for fit, sec
B0 = fit1.x0; % resonance center in Gauss

if ~exist('var0')
    var0 = [Ipp*0.005 B0 T2 Hpp*0.1 Ipp*0.005 B0 T2 Hpp*0.1];
end

%print guess parameters
var0

% function to minimize: sum of squared errors
fitfunc = @(var) abs(var(1))*ESRVoigtSimulation(x, var(2), 1e-20, abs(var(3)), Bmw, abs(var(4)), 1, Pars.B0MA*1e4)' + abs(var(5))*ESRVoigtSimulation(x, var(6), 1e-20, abs(var(7)), Bmw, abs(var(8)), 1, Pars.B0MA*1e4)';
sumofsquares = @(var) sum(sum( abs(fitfunc(var) - y).^2  ));

% Fit model to data with fminsearch (Nelder Mead algorithm, much better
% convergance than Levenberg Marquard or trust Region)
opt = optimset('TolFun',1e-12,'TolX',1e-12,'PlotFcns',@optimplotfval, 'MaxFunEvals', 1e6, 'MaxIter', 1e6);
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

    [Chi(i), dChi(i)] = SusceptebilityCalc(FittedArea(:,i), Pars);
    [NSpin(i), dNSpin(i)] = SpinCounting(FittedArea(:,i), Pars);
end

% and error margin
dChi = dArea .* Chi./FittedArea ;

% compute signal to noise
SNR = snr(yFit,yFit-y);

Bp2p = FWHMLorentz*1e4/sqrt(3);
dBp2p = dFWHMLorentz*1e4/sqrt(3);

% save data to output array
outputarray = {T, Pars.GFactor, B0 Brms, Bp2p, dBp2p, Chi, dChi, NSpin, dNSpin, SNR};
outputnames = {'T','gfactor', 'B0', 'Brms', 'Bp2p', 'dBp2p', 'Chi','dChi' ,'NSpin','dNSpin','SNR'};

for k=1:length(outputarray)
    output.(outputnames{k}) = outputarray(k);
end

output.x=x;
output.y=y;

output.T2 = T2;

close all

plot(x,y,'.');
hold on
plot(x, yFit);

end