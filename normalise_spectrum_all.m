function [dset] = normalise_spectrum_all(dset)
%NORMALISE_SPECTRUM_ALL normalises an EPR spectrum for all acquisition conditions
%
% 	ESR data is normalised for receiver gain, number of scans (Ns = 1),
% 	time constant (Tc = 1 ms). Additionally, we also normalise for Q-factor,
%   (Q = 10,000), modulation amplitude (ModAmp = 1 Gauss) and microwave
%   power (Pmw = 2 mW) so that we can compare spectra recorded under
%   different conditions.
%
%   Note: All of the additional parameters can influence both the signal
%         amplitude and shape. Use the results with caution.
%
%   INPUT:
%   dset - dataset from BrukerRead
%
% 	OUTPUT:
% 	dset - normalised dset

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$

import esr_analyses.*
import esr_analyses.utils.*

%% Normalise for default parameters first
[dset] = normalise_spectrum(dset);

%% Normalise for other parameters next
pars = dset.Properties.UserData;

QValue_norm = 10000;
B0MA_norm   = 1e-4;  % 1 Gauss
MWPW_norm   = 2e-3;  % 2 mW

if is_powersat_exp(pars) % don't normalise for MWPW in power saturation measurement
    dset{:, 2:end} = dset{:, 2:end} * QValue_norm/pars.QValue * B0MA_norm/pars.B0MA;
else
    dset{:, 2:end} = dset{:, 2:end} * QValue_norm/pars.QValue * B0MA_norm/pars.B0MA * sqrt(MWPW_norm)/sqrt(pars.MWPW);
end


end