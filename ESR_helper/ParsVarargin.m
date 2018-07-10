function [ argout ] = ParsVarargin( input )
%PARSVARARGIN Pairs up argument discriptions with values form varargin

if mod(length(input),2)~=0
    error('More parameters than arguments given');
end
argout = struct();
for i = 1:2:length(input)
    argout.(input{i}) = input{i+1};
end

end

