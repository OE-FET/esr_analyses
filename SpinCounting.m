function [NSpin, dNSpin, Data] = SpinCounting(varargin)
%SPINCOUNTING Performs spin-counting from epr spectrum
%   SYNTAX:
%   [NSpin,dNSpin,Data] = SPINCOUNTING(x, y, Pars)
%   [NSpin,dNSpin,Data] = SPINCOUNTING(Path2File)
%   [NSpin,dNSpin,Data] = SPINCOUNTING()
%   [NSpin,dNSpin,Data] = SPINCOUNTING(DoubleInt, Pars)
%
%   OUTPUT:
%   NSpin - number of spins
%   dNSpin - estimated error
%   Data.xNorm - x-axis data
%   Data.yCorr - baseline corrected ESR spectrum
%   Data.Int1 - integrated ESR spectrum
%   Data.Pars - experiment parameters
%
%   DEPENDENCIES:
%   NormaliseSpectrum.m
%   MWmean.m
%   DoubleIntNUM.m
%   Natural Constants
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2018/07/05 12:58 $    $Revision: 0.1 $

%% INPUT PROCESSING

narginchk(0, 4);

baseline = 'y';

switch nargin
    case 0
        [x, y, Pars] = NormaliseSpectrum;
    case 2
        DoubleInt    = varargin{1};
        Pars         = varargin{2};
    case 1
        [x, y, Pars] = NormaliseSpectrum(varargin{1});
    case 3
        [x, y, Pars] = NormaliseSpectrum(varargin{1}, varargin{2}, varargin{3});
    case 4
        [x, y, Pars] = NormaliseSpectrum(varargin{1},varargin{2},varargin{3});
        baseline     = varargin{4};
end

% check if file is a sliced 2D spectrum. If yes, use MW power from slice
if isfield(Pars, 'PROCESS') && strcmp(Pars.PROCESS, '''prSlice''')
    array = strsplit(Pars.Microwave);
    Pars.MWPW = str2double(array{end})*1e-3;
end

%% account for MW field distribution in cavity
try
    position_correction = MWmean(Pars);
catch
    fprintf(['MW field distribution in the cavity could not be read from the ' ...
        'DSC file.\nProceeding without correcting for the position of the '...
            'sample in the cavity.\n \n']);
    position_correction=1;
end

%% Calculate number of Spins

% Newer Bruker spectrometers save the cavity conversion factor in the DSC
% file. Try to read it out, otherwise prompt user for input. 
if isfield(Pars,'ConvFact')==0
    Pars.ConvFact = input('Please input the cavity specific conversion factor (from spin counting calibration):\n[default = 9.2710e-09]\n');
    if isempty(Pars.ConvFact)
        Pars.ConvFact = 9.2710e-09;
    end
end

% Get sample temperature from parameter file or promt user for input
if isfield(Pars, 'Temperature')==0
    T = input('Please give measurement temperature in K:\n[default = 298 K]\n');
    if isempty(T)
        T = 298;
    end
else
    T = str2double(strtrim(regexprep(Pars.Temperature,'K','')));
end

% get QValue from parameter file or promt user for input
if isfield(Pars, 'QValue')==0
    Pars.QValue=input('Please give cavity Q-value: ');
end

% integrate the ESR spectrum twice
if nargin ~=2
    [DoubleInt, SingleInt, yCorr] = DoubleIntNUM(x, y, baseline);
end

S=1/2;

boltzman_factor = planck*Pars.MWFQ / (boltzm*T);

% Get MW powers
if strcmp(Pars.YTYP, 'IGD')==1
    Pmw = Pars.z_axis/1000; % MW Power in W
else
    Pmw = Pars.MWPW;
end

% calculate the number of spins in the sample
% make sure to adjust the calibration factor k for every cavity

% Cavity and Bridge calibration factor
k = 200/(Pars.BridgeCalib * Pars.ConvFact);

NSpin = k * DoubleInt ./ (Pars.QValue * sqrt(Pmw) * Pars.B0MA * S*(S+1) * ...
    boltzman_factor * position_correction);

% estimate the error
dI = 0.04*DoubleInt;
dQ = 0.03*Pars.QValue;

dNSpin = dI.*NSpin./DoubleInt + dQ*NSpin/Pars.QValue;

% create ouput data structure
if nargin ~= 2
    Data.x          = x;
    Data.yCorr      = yCorr;
    Data.SingleInt  = SingleInt;
    Data.DoubleInt  = DoubleInt;
    Data.Pars       = Pars;
else
    Data = {};
end

end