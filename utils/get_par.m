function pars = get_par(pars, name, default)
%GET_PAR Fills out missing parameters by asking user for input.

if ~isstruct(pars); error('First argument must be a structure.'); end
if ~ischar(name); error('Second argument must be a string / character array.'); end

if ~isfield(pars, name) && nargin == 3
    
    msg = 'Please give value for parameter "%s" [default = %s]: ';
    
    pars.(name) = input(sprintf(msg, name, num2str(default)));
    
    if isempty(pars.(name))
        pars.(name) = default;
    end

elseif ~isfield(pars, name) && nargin == 2

    pars.(name) = [];
    
    while isempty(pars.(name))
        pars.(name) = input(['Please give value for parameter ' name ': ']);
    end
end

end