function [y, l] = multi_voigtian(Sys, Exp, Opt, varargin)
%MULTI_VOIGTIAN Multiple Pseudo-Voigt peaks with paraeters in Sys, Exp.
%
% 	MULTI_VOIGTIAN is designed to work in conjuction with the easyspin
%	toolbox. The structures Sys and Exp contain fitting parameters and
%	experimental conditions, respectively. The strutcure Opt is only a
%	dummy variable for compatability.

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2018/07/05 12:58 $    $Revision: 0.1 $

n = get_kwarg(varargin, 'deriv', 1);

x = linspace(Exp.Range(1), Exp.Range(2), Exp.nPoints);
l = zeros(Opt.n, length(x));

for n=1:Opt.n

    % calculate indivual voigt functions
    l(n,:) = Sys.Area(n) * pseudo_voigt(x, Sys.B0(n), ...
        Sys.FWHM_gauss(n), Sys.FWHM_lorentz(n), 'deriv', n);

end

% sum over all components
y = sum(l, 1);

end