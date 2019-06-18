function argout = plot2D(x, y, o, varargin)

style    = get_kwarg(varargin, 'style', '');

[X, Y] = meshgrid(x, y);

if ~all(size(o) == [length(y) length(x)])
    o = reshape(o, [length(x) length(y)])';
end

if ~all(size(o) == [length(y) length(x)])
    error('Dimensions of axes and ordinate data do not match.')
end

h = contourf(X, Y, o, style);

colorbar;

if nargout > 1
    argout = h;
end

end