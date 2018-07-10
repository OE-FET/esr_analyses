function [xNorm,yNorm,Pars]=NormaliseSpectrumAll(varargin)
%NORMALISESPECTRUMALL normlizes an EPR spectrum for all aquisition conditions
%
% 	ESR data is normalised for reciever gain, number of scans (Ns = 1),
% 	time constant (Tc = 1 ms). Aditionally, we also normalise for Q-factor,
% 	modulation amplitude and microwave power so that we can compare spectra
% 	recorded under different conditions.
%
% 	INPUT:
% 	NORMALISESPECTRUMALL() promts user for path input via GUI
% 	NORMALISESPECTRUMALL(Path) loads data from specified path
% 	NORMALISESPECTRUMALL(x,y,Pars) uses given data in (x,y) and experimental
% 	conditions from Pars
%
% 	OUTPUT: 
% 	x_norm and y_norm are vectors containg the magnetic field in Gauss and the
% 	normalised ESR signal intensity.
%
% 	DEPENDENCIES:
% 	NormaliseSpectrum.m
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2018/07/05 12:58 $    $Revision: 0.1 $

%% Input analyses
% load file
narginchk(0,3);
switch nargin
    case 0
        [x,y,Pars] = NormaliseSpectrum;
    case 1
        [x,y,Pars] = NormaliseSpectrum(varargin{1});
    case 3
        x    = varargin{1};
        y    = varargin{2};
        Pars = varargin{3};
        [x,y,Pars] = NormaliseSpectrum(x,y,Pars);
end

% Does the file contain multiple scans at different powers (2D) or only one scan?
if strcmp(Pars.YTYP, 'IGD')==1
    % don't normalise for MWPW in power saturation measurement
    yNorm = y / (Pars.QValue * Pars.B0MA );
else
    yNorm = y / (Pars.QValue * sqrt(Pars.MWPW) * Pars.B0MA );
end

%% normalise y-axis
xNorm=x;

end