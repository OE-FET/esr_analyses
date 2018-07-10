function [output] = ESRAnalysesNUM(varargin)
%ESRANALYESNUM performs normalization and spin-counting of an ESR signal by
%numercial double integration.
%   The ESR signal is normalised according to measurement conditions. If
%   required, a background signal can be subtracted before performing the
%   analyses.
%
%   The total ESR intensity and spin susceptibility are determined by
%   numerical double-integration. The number of spins is then calculated 
%   by assuming that the sample follows a Curie-Weiss law. This assumption is
%   valid when T >> E/2. X-Band EPR typically operates around B = 350 mT while
%   a field of up to 0.5 T corresponds to E/2 = 28.9 ueV. This still remains
%   far lower than the corresponding thermal energy of 450 ueV at 5 K. 
%
%   INPUT(S):
%   ESRAnalysesNUM()            - prompts user for spectrum file
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
%
%   DEPENDENCIES:
%   SubtractBackground.m
%   NormaliseSpectrum.m
%   MarkerCalib.m
%   gfactor_determination.m
%   SpinCounting.m
%   num2clip.m
%   getSamplePosition.m
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2018/07/05 12:58 $    $Revision: 1.1 $

% load and normalise spectrum, subtract a background spectrum if requested
close all

switch nargin
    case 3
        x = varargin{1}; y = varargin{2}; Pars = varargin{3};
    case  2
        [x,y,Pars]=SubtractBackground(varargin{1},  varargin{2});
    case 1
        str = input('Would you like to subtract a background signal? y/[n]','s');
        if strcmp(str,'y')
            [x,y,Pars]=SubtractBackground(varargin{1});
        else
            [x,y,Pars]=NormaliseSpectrum(varargin{1});
        end
    case 0
        str = input('Would you like to subtract a background signal? y/[n]','s');
        if strcmp(str,'y')
        [x,y,Pars]=SubtractBackground;
        else
        [x,y,Pars]=NormaliseSpectrum;
        end
end
        

% if a marker is used for g-factor calibration, normalise x-axis according
% to marker position
% this function does nothing if no marker signal is detected
[x,y,Pars]=MarkerCalib(x,y,Pars);
sample_g=gfactor_determination(x,y,Pars);

% try to load sample temperature from parameter file
% prompt user for entry if not found
try
    T = str2double(strtrim(regexprep(Pars.Temperature,'K','')));
catch
    T = input('Please give the measurement temperature in Kelvin:');
    Pars.Temperature = [num2str(T), ' K'];
end

% count the number of spins, accounting for MW field distribution in cavity
[NSpin, dNSpin, Data] = SpinCounting(x,y,Pars);

% calculate suscepebility
S = 1/2; mu0 = 4*pi*10^(-7);
Chi = NSpin .* (mu0 * S*(S+1) * sample_g.^2 * bmagn^2)/(3*boltzm*T );

% and error margin
dChi = dNSpin*Chi/NSpin;

% save data to output array
outputarray={T, 100/T, sample_g, Chi, dChi, NSpin, dNSpin};
outputnames={'T','InverseT','g', 'Chi','dChi' ,'NSpin','dNSpin'};

for k=1:length(outputarray)
    output.(outputnames{k})=outputarray(k);
end

% save output array to base workspace
assignin('base', 'outputarray', outputarray);
% save output array to clipboard
copy(outputarray);

output.xNorm=x;
output.yNorm=Data.yCorr;

end