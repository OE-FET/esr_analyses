function [x_norm, y_norm, pars] = normalise_spectrum_all(x, y, pars)
%NORMALISE_SPECTRUM_ALL normlizes an EPR spectrum for all aquisition conditions
%
% 	ESR data is normalised for reciever gain, number of scans (Ns = 1),
% 	time constant (Tc = 1 ms). Aditionally, we also normalise for Q-factor,
% 	modulation amplitude and microwave power so that we can compare spectra
% 	recorded under different conditions.
%
% 	INPUT:
% 	normalise_spectrum_all(x, y, pars) uses given data in (x,y) and 
% 	experimental conditions from pars.
%
% 	OUTPUT: 
% 	x_norm and y_norm are vectors containg the magnetic field in Gauss and
% 	the normalised ESR signal intensity.
%
% 	DEPENDENCIES:
% 	normalise_spectrum.m
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2018/07/05 12:58 $    $Revision: 0.1 $

%% Normalise for default parameters first
[x, y, pars] = normalise_spectrum(x, y, pars);

%% Normalise for other parameters next

% Does the file contain multiple scans at different powers (2D) or only one scan?
if strcmp(pars.YTYP, 'IGD') == 1
    % don't normalise for MWPW in power saturation measurement
    y_norm = y / (pars.QValue * pars.B0MA );
else
    y_norm = y / (pars.QValue * sqrt(pars.MWPW) * pars.B0MA );
end

%% normalise y-axis
x_norm=x;

end