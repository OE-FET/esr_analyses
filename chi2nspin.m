function [n_spin] = chi2nspin(chi, T, S)
%CHI2NSPIN Convert the Curie susceptebility to a number of spins.
%
%   SYNTAX:
%   [n_spin] = CHI2NSPIN(chi)
%   [n_spin] = CHI2NSPIN(chi, T)
%   [n_spin] = CHI2NSPIN(chi, T, S)
%
%   INPUT(S):
%   chi - magnetc susceptibility in SI units
%   T - temperature in Kelvin (default: T = 298)
%   S - electron spin (default: S = 1/2)
%
% 	OUTPUT(S):
%   n_spin - number of spins per cubic meters
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2019/05/06 12:58 $    $Revision: 1.1 $

%% input analyses
if nargin < 3; S=1/2; end
if nargin < 2; T=298; end

%% Susceptebility calculation
n_spin =  chi * ( mu0 * S*(S+1)*gfree^2*bmagn^2 / (3 * boltzm *T)  )^(-1);

end