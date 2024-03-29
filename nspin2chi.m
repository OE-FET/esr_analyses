function [chi] = nspin2chi(n_spin, T, S)
%CHI2NSPIN Convert the Curie susceptibility to a number of spins.
%
%   SYNTAX:
%   [chi] = CHI2NSPIN(n_spin)
%   [chi] = CHI2NSPIN(n_spin, T)
%   [chi] = CHI2NSPIN(n_spin, T, S)
%
%   INPUT:
%   n_spin - number of spins per cubic meters
%   T - temperature in Kelvin (default: T = 298)
%   S - system spin (default: S = 1/2)
%
% 	OUTPUT:
%   chi - magnetc susceptibility in SI units (i.e., dimensionless)
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$

import esr_analyses.*
import esr_analyses.utils.*

%% input analyses
if nargin < 3; S=1/2; end
if nargin < 2; T=298; end

%% Susceptebility calculation
chi = n_spin * mu0 * S*(S+1)*gfree^2*bmagn^2 / (3*boltzm *T);

end

