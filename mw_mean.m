function [position_correction] = mw_mean(pars)
%MW_MEAN Accounts for non-uniform MW field distribution in cavity by
% calculating effective "f_mean" depending on the vertical position and
% length of sample in cavity. The field distribution is calculated from the
% polynomian coefficients saved in the .DSC file by Xepr.
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$

import esr_analyses.*
import esr_analyses.utils.*

% PolyCof contains the polynomial coefficients from lowest to
% higherst order. However, Matlab convention goes from highest to lowest.
% we need to flip PolyCof.
coefVector = fliplr(pars.PolyCof(2,:));

% read out cavity dimensions
ResCenter = str2double(strtrim(regexprep(pars.ResCenter, 'mm', '')));
ResLength = str2double(strtrim(regexprep(pars.ResLength, 'mm', '')));

% ask for height and length of sample if not in Pars
pars = get_sample_position(pars);

% check if sample length is smaller than the cavity size
% if length > ResLength, prompt user to enter new value 3 times, then abort
count = 0;
while pars.SampleL > ResLength
    disp('Sample length is larger than the cavity height. Please enter a valid sample length.');
    pars.SampleL = input('Please give sample length in mm:\n[default = 20 mm]\n');
    if isempty(pars.SampleL); pars.SampleL=25; end
    count = count+1;
    if count > 3; error('Invalid sample length.'); end
end

% Calculate average MW magnetic field over sample in cavity
sampleDim = ((pars.SampleH-pars.SampleL/2):0.01:(pars.SampleH+pars.SampleL/2));

B1 = polyval(coefVector, sampleDim - ResCenter);
position_correction = mean(B1);


end

