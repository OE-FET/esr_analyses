function [argout, input] = get_varargin(input, name, default)
%GET_VARARGIN Pairs up argument discriptions with values form varargin.

for i = 1:length(input)
    if ischar(input{i}) && strcmp(input{i}, name)
        argout = input{i+1};
        input(:,i+1) = [];
        input(:,i) = [];
        return
    end
end

if ~exist('argout', 'var') && nargin == 3
    argout = default;
end

end

