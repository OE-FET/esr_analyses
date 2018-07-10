function vp = Voigt_fadf(x, gam, sig)

% VOIGT returns Voigt profile data.
% x             Vector of x values.
% x0            Center of voigt function
% gammaVal      Gamma value.
% FWHMGauss          Sigma value.


% converting to dimensionless coordinates

z = arrayfun(@(q) (q+1i*gam)/(sqrt(2)*sig), x);
vp = (1/(sig*sqrt(2*pi))) * real(fadf(z)); % Get Voigt from Faddeeva fn.
vp = vp./max(vp);


end