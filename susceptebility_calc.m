function [Chi, dChi] = susceptebility_calc(DoubleIntArea, Pars)
%SUSCEPTEBILITY_CALC calculates the magnetic susceptebility from given 
% integrated area of an EPR spectrum and measurement parameters
%
%   SYNTAX:
%   [Chi, dChi, Data] = susceptebility_calc(DoubleIntArea, Pars)
%
%   INPUT(S):
%   DoubleIntArea - Integrated area of absorption EPR spectrum or double 
%       integral from 1st harmonic spectrum
%   Pars - structure containing measurement parameters from EPR experiment
%
%   OUTPUT(S):
%   Chi - magnetic susceptebility
%   dChi - estimated error
%
%   DEPENDENCIES:
%   Natural Constants
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2018/07/05 12:58 $    $Revision: 0.1 $

%% INPUT PROCESSING

% check if file is a sliced 2D spectrum. If yes, use MW power from slice
if isfield(Pars, 'PROCESS') && strcmp(Pars.PROCESS, '''prSlice''')
    array = strsplit(Pars.Microwave);
    Pars.MWPW = str2double(array{end})*1e-3;
end

%% account for MW field distribution in cavity
try
    position_correction = mw_mean(Pars);
catch
    fprintf(['MW field distribution in the cavity could not be read from the ' ...
        'DSC file.\nProceeding without correcting for the position of the '...
            'sample in the cavity.\n \n']);
    position_correction=1;
end

%% Calculate susceptebility

% Newer Bruker spectrometers save the cavity conversion factor in the DSC
% file. Try to read it out, otherwise prompt user for input. 
if isfield(Pars, 'ConvFact')==0
    Pars.ConvFact = input('Please input the cavity specific conversion factor (from spin counting calibration):\n[default = 9.2710e-09]\n');
    if isempty(Pars.ConvFact)
        Pars.ConvFact = 9.2710e-09;
    end
end

% get QValue from parameter file or promt user for input
if isfield(Pars, 'QValue')==0
    Pars.QValue = input('Please give cavity Q-value: ');
end

% get g-factor from parametr file
sample_g = Pars.GFactor;

% get MW powers
if strcmp(Pars.YTYP, 'IGD')==1
    Pmw = Pars.z_axis/1000; % MW Power in W
else
    Pmw = Pars.MWPW;
end

% make sure to adjust the calibration factor k for every cavity

% get cavity and Bridge calibration factor
k = 200/(Pars.BridgeCalib * Pars.ConvFact);

% normalization factor for measurement conditions
norm = (Pars.QValue * sqrt(Pmw) * Pars.B0MA * position_correction);

Chi = mu0* k * DoubleIntArea .* sample_g.^2 .* bmagn^2 ./ (3 * planck * Pars.MWFQ  * norm);

% estimate the error
dI   = 0.04*DoubleIntArea;
dQ   = 0.03*Pars.QValue;
dChi = dI.*Chi./DoubleIntArea + dQ.*Chi./Pars.QValue;

end
