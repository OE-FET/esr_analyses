function [dset] = normalise_spectrum_all(dset)
%NORMALISE_SPECTRUM_ALL normlizes an EPR spectrum for all aquisition conditions
%
% 	ESR data is normalised for reciever gain, number of scans (Ns = 1),
% 	time constant (Tc = 1 ms). Aditionally, we also normalise for Q-factor,
%   (Q = 10,000), modulation amplitude (ModAmp = 1 Gauss) and microwave
%   power (Pmw = 1mW) so that we can compare spectra recorded under
%   different conditions.
%
%   Note: All of the additional parameters can influence both the signal
%         amplitude and shape. Use the results with caution.
%
%   INPUT(S):
%   dset from BrukerRead
%
% 	OUTPUT(S):
% 	normalised dset

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$

import esr_analyses.*
import esr_analyses.utils.*

%% Normalise for default parameters first
[dset] = normalise_spectrum(dset);

%% Normalise for other parameters next
pars_norm = pars;
pars_norm.QValue = 1;
pars_norm.B0MA = 1;

dset.Properties.UserData.QValue = 10000;
dset.Properties.UserData.B0MA = 1;

if is_powersat_exp(pars) % don't normalise for MWPW in power saturation measurement
    dset{:, 2:end} = dset{:, 2:end} / (pars.QValue * pars.B0MA);
else
    pars_norm.MWPW = 1;
    dset{:, 2:end} = dset{:, 2:end} / (pars.QValue * sqrt(pars.MWPW) * pars.B0MA );
end


end