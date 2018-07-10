function nDiff = comparePars(Pars1, Pars2)
%COMPAREPARS Compares experimental conditions ffom two ESR measurements
%   nDiff = COMPAREPARS(Pars1, Pars2) compares the measurement parameters
%   Pars1 and Pars2 from two different ESR data sets and outputs the number
%   of different parameters.
%
		
%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2018/07/05 12:58 $    $Revision: 1.1 $

%% determine parameters to be ignored in comparision
% in case of a power dependant measurement, ignore MWPW paramaters

if strcmp(Pars1.YTYP, 'IGD') == 1
    ParCheckIgnore = {'TITL', 'TIME', 'MWFQ', 'FrequencyMon', 'Flyback', 'Power', 'PowerAtten', 'MWPW'};
else
    ParCheckIgnore = {'TITL', 'TIME', 'MWFQ', 'FrequencyMon', 'Flyback'};
end

ParsNames1 = fieldnames(Pars1);
ParsNames2 = fieldnames(Pars2);

N = min(length(ParsNames1), length(ParsNames2));

comp = zeros(N, 1);

for i = 1:N
    if isequal(ParsNames1(i), ParsNames2(i)) && ~StrInList(ParsNames1(i), ParCheckIgnore)
        comp(i) = isequal(Pars1.(ParsNames1{i}), Pars2.(ParsNames2{i})) - 1;
    end
end

nDiff = find(comp);
if nDiff>0
    fprintf(2, 'The following parameters are different:\n');
end
for i = nDiff'
    disp([ParsNames1(i), Pars1.(ParsNames1{i}), 'vs', Pars2.(ParsNames2{i})]);
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

