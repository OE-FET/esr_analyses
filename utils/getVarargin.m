function [argout, input] = getVarargin(input, nameStr)
%PARSVARARGIN Pairs up argument discriptions with values form varargin

for i = 1:length(input)
    if ischar(input{i}) && strcmp(input{i}, nameStr)
        argout = input{i+1};
        input(:,i+1) = [];
        input(:,i) = [];
        return
    end
end

end

