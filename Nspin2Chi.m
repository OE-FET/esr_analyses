function [chi] = nspin2chi(n_spin, T, S)
%CHI2NSPIN Convert the Curie susceptebility to a number of spins.
%   [chi] = nspin2chi(n_spin)
%   [chi] = nspin2chi(n_spin, T)
%   [chi] = nspin2chi(n_spin, T, S)
%
%   INPUT(S):
%   n_spin - number of spins per cubic meters
%   T - temperature in Kelvin (default: T = 298)
%   S - electron spin (default: S = 1/2)
%
% 	OUTPUT(S):
%   chi - magnetc susceptibility in SI units
%
% 	Dependencies:
% 	natural constants
%

%% input analyses
if nargin < 3; S=1/2; end
if nargin < 2; T=298; end

%% Susceptebility calculation
chi = n_spin * mu0 * S*(S+1)*gfree^2*bmagn^2 / (3*boltzm *T);

end

