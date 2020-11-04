function n_diff = compare_pars(dset1, dset2)
%COMPARE_PARS Compares experimental conditions from two ESR measurements
%
%   n_diff = COMPARE_PARS(dset1,dset2) compares the measurement parameters
%   from two different ESR data sets dset 1 and dset2 and outputs the
%   number of different parameters. It also prints the different parameter
%   values to the console.
%
%   Instead of dataset tables, this function can also take parameter
%   structures as input.
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$

%% determine parameters to be ignored in comparision
% in case of a power-dependent measurement, ignore MWPW paramaters

import esr_analyses.*
import esr_analyses.utils.*

if istable(dset1)
    dset1 = dset1.Properties.UserData;
end
if istable(dset2)
    dset2 = dset2.Properties.UserData;
end

ParCheckIgnore = {'TITL', 'TIME', 'DATE', 'MWFQ', 'FrequencyMon', ...
         'Flyback', 'XMIN', 'XWID', 'XMAX', 'StaticFieldMon', ...
         'NbScansAcc', 'NbScansToDo', 'OPER', 'Stability'};

if strcmp(dset1.YTYP, 'IGD') == 1
    ParCheckIgnore = [ParCheckIgnore, {'Power', 'PowerAtten', 'MWPW'}];
end

ParsNames1 = fieldnames(dset1);
ParsNames2 = fieldnames(dset2);

ParsNames = intersect(ParsNames1, ParsNames2);
N = length(ParsNames);

comp = zeros(N, 1);

for i = 1:N
    if ~any(strcmp(ParCheckIgnore, ParsNames{i}))
        comp(i) = ~isequal(dset1.(ParsNames{i}), dset2.(ParsNames{i}));
    end
end

% print all parameters that are different
diff_indices = find(comp);
n_diff = length(diff_indices);
if n_diff > 0
    fprintf(2, 'The following parameters are different:\n');
end
for i = diff_indices'
    fprintf('%s:\t%s vs %s\n', ParsNames{i}, dset1.(ParsNames{i}), dset2.(ParsNames{i}));
end

end
