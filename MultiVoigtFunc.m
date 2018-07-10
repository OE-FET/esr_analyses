function [y, l] = MultiVoigtFunc(Sys, Exp, Opt)
%MULTIVOIGTFUNC Three Pseudo-Voigt peaks with paraeters in Sys, Exp
%
% 	MULTIVOIGTFUNC is designed to work in conjuction with the easyspin toolbox.
%	The structures Sys and Exp contain fitting parameters and experimental
%	conditions respectively. The strutcure Opt is only a dummy variable for
%	compatability.

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2018/07/05 12:58 $    $Revision: 0.1 $

x = linspace(Exp.Range(1), Exp.Range(2), Exp.nPoints);

Ap2p = zeros(Opt.n, 1);
l = zeros(Opt.n, length(x));

for n=1:Opt.n
    
    % convert areas to peak2peak amlitudes
    Ap2p(n) = -24*Sys.Area(n)/((8*sqrt(3)*pi*(Sys.s(n)-1) - 3*sqrt(2*exp(1)*pi)*Sys.s(n))*Sys.lwpp(n)^2);
    % calculate indivual Voigt functions
    l(n, :) = PseudoVoigtDeriv(x, Ap2p(n), Sys.lwpp(n), Sys.B0(n), Sys.s(n));
    
end
% sum over all components
y = sum(l,1);

end