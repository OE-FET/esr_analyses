function Bmw = get_mw_fields(pars)
%GET_MW_FIELDS Gets microwave fields from measurement parameter structure.
%
%   Converts MW powers given in pars structure to MW magentic field
%   amplitudes, given the measurement conditions from pars. In case of
%   power saturation measurements, a vector is returned. Otherwise,
%   GET_MW_FIELDS returns a scalar.
%
%   INPUT(S):
%   pars - structure containing measurement conditions, sample dimensions,
%          etc.
%
%   OUTPUT(S):
%   Bmw  - Microwave magnetic field amplitude(s) in Tesla
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2019/05/06 12:58 $    $Revision: 1.1 $

%% Get microwave powers based on whether this is a PowerSat measurement
if isfield(pars,'YNAM') & strcmp(pars.YNAM, '''Microwave Power''')
    Pmw = pars.z_axis*1e-3; % MW Power in W
else
    Pmw = pars.MWPW; % MW Power in W
end

% convert MWPW to magnetic field strength in Tesla
% use Qref = 7500 and c = 2.0 for Nagoya files
Qref    = 8355;
c       = 2.2;      % without cryostat mounted
cCryo   = 2.2*1.3;  % with crystat, calibrated with gMarker

pars    = get_sample_position(pars);
f_mean  = mw_mean(pars);
Bmw     = f_mean * cCryo * sqrt(Pmw) * sqrt(pars.QValue/Qref) * 1E-4; % in Tesla

end