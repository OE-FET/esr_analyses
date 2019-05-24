function [returnValue, output] = get_varargin(input, name, defaultValue)
%GET_VARARGIN Pairs up argument discriptions with values form varargin.

output = input;

for i = 1:length(output)
    if ischar(output{i}) && strcmp(output{i}, name)
        returnValue = output{i+1};
        output(:,i+1) = [];
        output(:,i) = [];
        return
    end
end

if ~exist('argout', 'var') && nargin == 3
    returnValue = defaultValue;
end

end

