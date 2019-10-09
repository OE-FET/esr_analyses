function pars = get_sample_position(pars)
%GET_SAMPLE_POSITION gets the sample position :)
%
% 	Ask the user for the length of the sample and the position of the
% 	sample in the ESR cavity (if they are not already present in the pars
%   structure. These are necessary for proper determination of the mw field
%   over the sample.
%
%   INPUT(S):
%   get_sample_position(pars)   - experimental parameters
%
% 	OUTPUT(S):
% 	pars.SampleL  - length of the sample [mm]
% 	pars.SampleH  - height of the sample inside the cavity [mm]
%
%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$

import esr_analyses.*
import esr_analyses.utils.*

if ~isfield(pars,'SampleL')

    pars.SampleL = input('Please give sample length in mm [default = 20]: ');
    if isempty(pars.SampleL)
        pars.SampleL = 20;
    end
end

if ~isfield(pars,'SampleH')

    pars.SampleH = input('Please give vertical position of sample center in the cavity,\nmeasured from the top collar in mm [default = 62.5]: ');
    if isempty(pars.SampleH)
        pars.SampleH = 62.5;
    end
end

end