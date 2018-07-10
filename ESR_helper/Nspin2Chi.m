function [ Chi ] = Nspin2Chi( varargin)
%NSPIN2CHI Convert number of spins to Curie susceptebility
%   [ Chi ] = Nspin2Chi( Nspin, T, S )
% T and S are optional arguments. If not given, the user is promted to
% input the temperature and S = 1/2 is used.
%
% Dependencies:
% natural constants package

%% input analyses
if nargin <1
    error('Please give the number of spins as input argument.')
elseif nargin < 2
    Nspin = varargin{1};
    T = input('Please give the sample temperature in K:');
    S = 1/2;
elseif nargin < 3
    Nspin = varargin{1};
    T = varargin{2};
    S = 1/2;
elseif nargin < 4
    Nspin = varargin{1};
    T = varargin{2};
    S = varargin{3};
else
    error('Too many arguments given.')
end

%% Susceptebility calculation
Chi = Nspin * mu0 * S*(S+1)*gfree^2*bmagn^2 / (3*boltzm *T);

end

