function [position_correction] = mw_mean(Pars)
%MW_MEAN Accounts for non-uniform MW field distribution in cavity by
% calculating effective "f_mean" depending on the vertical position and
% length of sample in cavity. The field distribution is calculated from the
% polynomian coefficients saved in the .DSC file by Xepr.
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2019/05/06 12:58 $    $Revision: 1.1 $

import esr_analyses.*
import esr_analyses.utils.*

%%
% read string with polymer coefficients for Bmw distribution
[info, coef] = strtok(Pars.PolyCof);
% coef contains the coefficient vector
coef = str2double(strsplit(strtrim(coef), ','));
% info contains information about how to read coef
info = str2double(strsplit(regexprep(info, {'{', '}'}, ''), {',', ';'}));

% reshape coef vector into matrix
coefMatrix = reshape(coef, info(2:3));
coefVector = coefMatrix(info(1), :);

% coefVector now contains the polynomial coefficients from lowest to
% higherst order. However, Matlab convention goes from highest to lowest.
% we need to flip coefVector.
coefVector=fliplr(coefVector);

% read out cavity dimensions
ResCenter = str2double(strtrim(regexprep(Pars.ResCenter, 'mm', '')));
ResLength = str2double(strtrim(regexprep(Pars.ResLength, 'mm', '')));

% ask for height and length of sample if not in Pars
Pars = get_sample_position(Pars);

% check if sample length is smaller than the cavity size
% if length > ResLength, prompt user to enter new value 3 times, then abort
count = 0;
while Pars.SampleL > ResLength
    disp('Sample length is larger than the cavity height. Please enter a valid sample length.');
    Pars.SampleL = input('Please give sample length in mm:\n[default = 20 mm]\n');
    if isempty(Pars.SampleL);Pars.SampleL=25;end
    count = count+1;
    if count > 3; error('Invalid sample length.'); end
end

% Calculate average MW magnetic field over sample in cavity
sampleDim = ((Pars.SampleH-Pars.SampleL/2):0.01:(Pars.SampleH+Pars.SampleL/2));

B1 = polyval(coefVector, sampleDim - ResCenter);
position_correction = mean(B1);


end

