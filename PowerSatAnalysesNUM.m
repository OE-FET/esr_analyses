function [argout, DATA, fitresult] = PowerSatAnalysesNum(varargin)
%POWERSATANALYSESNUM Numercial analyses of ESR power saturation curves
%
% 	Calculates power stauration curves from integrated intensities and maximum
% 	intensities. This programm numerically intergrates the ESR signals.
%
% 	Disadvantage: Long tails of the resonance peak are negletcted.
% 	Advantage: Works with every peak shape.
%
% 	Power saturation curves are given in DATA with the colmuns:
% 	[Pmw MwB II I_max DeltaBpp g]
%
% 	All spectra are collected in the structure 'argout.ERSIntensity'.
%
%	INPUTS:
%	POWERSATANALYSESNUM()							- pens GUI for file selection
%	POWERSATANALYSESNUM(x,y,Pars)					- data given by (x,y,Pars)
%	...('/path/to/signal')							- reads data from file
%	...('/path/to/signal', '/path/to/background')	- reads data and background from file
%
%	OUTPUT:
%	argout  - structure containing all fitting results and measurement data
%	DATA 	- matrix with power saturation curves
%	fitresult - structure containing fit results
%
% 	DEPENDENCIES:
% 	BrukerRead.m
% 	normalise_spectrum.m
% 	subtract_background.m
% 	gfactor_determination.m
% 	mw_mean.m
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2018/07/05 12:58 $    $Revision: 0.1 $

%% Read data 
if nargin == 1
    Path = varargin{1};
    str = input('Would you like to subtract a background signal? y/[n]', 's');
    if strcmp(str, 'y') == 1
        [x, y, pars] = subtract_background(Path);
    else
        [x, y, pars] = BrukerRead(Path);
    end
elseif nargin == 0
    str = input('Would you like to subtract a background signal? y/[n]','s');
    if strcmp(str, 'y') == 1
    [x, y, pars] = subtract_background;
    else
    [x, y, pars] = BrukerRead();
    end
elseif nargin == 3
    x = varargin{1}; y = varargin{2}; pars = varargin{3};
end

% Get MW powers
if strcmp(pars.YTYP, 'IGD') == 1
    Pmw = pars.z_axis/1000; % MW Power in mW
else
    error('The specified file is not a 2D data file.');
end

% get Q value
try
    QValue = pars.QValue;
catch
    QValue = input('Please give cavity QValue:');
end

%% calculate integrated intensities
% normalise for measurement conditions
[x, y] = normalise_spectrum(x, y, pars);

% perform double integration of all spectra
baseline = input('Perform base-line correction individually or as batch? i/[b]?','s');
if baseline == 'i'
    for i = 1:size(y, 2)
        [DoubleInt(i), ~, yCorr(:,i)] = double_int_num(x, y(:,i), 'y');
    end
    DoubleInt = DoubleInt';
else
    [DoubleInt, ~, yCorr] = double_int_num(x, y, 'y');
end

% detremine g-factors of all spectra
g=gfactor_determination(x, y(:,round(end/2)), pars);

% convert MWPW to magnetic field strength in Tesla
% use Qref = 7500 and c = 2.0 for Nagoya files!!!
% appropriate conversion factors and Qref values are availavbe form Bruker

Qref = 8355;
conv = 2.2;
cCryo = 2.2*1.3; % calbirated with gMarker

f_mean = mw_mean(pars);

Bmw = f_mean * cCryo * sqrt(Pmw) * sqrt(QValue/Qref) * 10^(-4); % in Tesla

% collect values in output matrix for easy copying to Origin
DATA = [Pmw Bmw DoubleInt];

xNorm = x;
yNorm = yCorr;

argout.MircowavePower = DATA(:,1);
argout.MicrowaveB = DATA(:,2);
argout.DoubleIntegrated = DATA(:,3);
argout.gValue = g;
argout.xNorm = xNorm;
argout.yNorm = yNorm;

%% Fit and plot results

try
fitresult = saturation_area_fit(argout.MicrowaveB(1:end), ...
    abs(argout.DoubleIntegrated(1:end)), argout.gValue);
catch
    disp('Fit could not be performed.');
end