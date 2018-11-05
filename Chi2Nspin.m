function [ Nspin ] = Chi2Nspin( varargin)
%NSPIN2CHI Convert the Curie susceptebility to anumber of spins. 
%   [ Chi ] = Chi2Nspin( Nspin, T, S )
% 	T and S are optional arguments. If not given, the user is promted to
% 	input the temperature and S = 1/2 is used.
%
% 	Dependencies:
% 	natural constants

%% input analyses
if nargin <1
    error('Please give the number of spins as input argument.')
elseif nargin < 2
    Chi = varargin{1};
    T = input('Please give the sample temperature in K:');
    S = 1/2;
elseif nargin < 3
    Chi = varargin{1};
    T = varargin{2};
    S = 1/2;
elseif nargin < 4
    Chi = varargin{1};
    T = varargin{2};
    S = varargin{3};
else 
    error('Too many arguments given.')
end

%% Susceptebility calculation
Nspin =  Chi * ( mu0 * S*(S+1)*gfree^2*bmagn^2 / (3 * boltzm *T)  )^(-1);

end
