function n_diff = compare_pars(pars1, pars2)
%COMPARE_PARS Compares experimental conditions from two ESR measurements
%
%   nDiff = COMPARE_PARS(pars1,pPars2) compares the measurement parameters
%   Pars1 and Pars2 from two different ESR data sets and outputs the number
%   of different parameters.
%
		
%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2019/05/06 12:58 $    $Revision: 1.1 $

%% determine parameters to be ignored in comparision
% in case of a power-dependent measurement, ignore MWPW paramaters

ParCheckIgnore = {'TITL', 'TIME', 'DATE', 'MWFQ', 'FrequencyMon', ...
         'Flyback', 'XMIN', 'XWID', 'XMAX', 'StaticFieldMon', ...
         'NbScansAcc', 'NbScansToDo', 'OPER'};

if strcmp(pars1.YTYP, 'IGD') == 1
    ParCheckIgnore = [ParCheckIgnore, {'Power', 'PowerAtten', 'MWPW'}];
end

ParsNames1 = fieldnames(pars1);
ParsNames2 = fieldnames(pars2);

N = min(length(ParsNames1), length(ParsNames2));

comp = zeros(N, 1);

for i = 1:N
    if isequal(ParsNames1(i), ParsNames2(i)) && ~StrInList(ParsNames1(i), ParCheckIgnore)
        comp(i) = isequal(pars1.(ParsNames1{i}), pars2.(ParsNames2{i})) - 1;
    end
end

% print all parameters that are different
n_diff = find(comp);
if n_diff>0
    fprintf(2, 'The following parameters are different:\n');
end
for i = n_diff'
    disp([ParsNames1(i), pars1.(ParsNames1{i}), 'vs', pars2.(ParsNames2{i})]);
end

end


function out = StrInList(str, list)
out = 0;
for j = 1:length(list)
    if isequal(str, list(j))
       out = 1;
    end
end
end

