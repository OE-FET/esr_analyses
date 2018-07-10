function [ xNorm, Pars ] = normMWFQ( x,Pars )

% Normalise magnetic field axis to a MWFQ of 9.6 GHz
offsetB = 10000*(9.6e9-Pars.MWFQ)*planck/bmagn/2;
xNorm = x+offsetB;

% update parameter set accordingly
Pars.MWFQ = 9.6*10^9;

end

