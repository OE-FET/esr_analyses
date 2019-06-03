function [returnValue, output] = get_kwarg(input, name, defaultValue)
%GET_KWARG Pairs up argument discriptions with values form varargin.

output = input;

% get parameter from varargin 'input' and return varargin without patameter
% as 'output'
for i = 1:length(output)
    if ischar(output{i}) && strcmp(output{i}, name)
        returnValue = output{i+1};
        output(:,i+1) = [];
        output(:,i) = [];
        return
    end
end

% if parameter could not be found, use default or prompt user for input
if ~exist('returnValue', 'var')
    if nargin == 3
        returnValue = defaultValue;
    else
        returnValue = [];
        while isempty(returnValue)
            msg = sprintf('Please give input parameter "%s": ', name);
            returnValue = input(msg);
        end
    end
end

end
