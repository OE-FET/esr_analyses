function [IArea, SingleInt, yFit, FitParam,FitError] = DoubleIntFIT(x,y)
%DOUBLEINTFIT Fits spectrum with a Voigt function and then performs double 
%integration on the fit
%   [IArea, SingleInt, yFit, FitParam, FitError] = DOUBLEINTFIT(x, y) fits 
%   an ESR spectrum with a Voigt function and then performs a double
%	integration on the fitted curve.
%
%   INPUT:
%   (x,y) - spectrum data
%
%   OUTPUT:
%   IArea - total double integrated area
%   Int1 - curve after first integration
%   fitresult - fitted Voigtian curve
%   fiterrors - errors from Voigtian fit
%
%   DEPENDENCIES:
%   VoigtFit.m
%   Voigt.m
% 

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2018/07/05 12:58 $    $Revision: 1.1 $

% determine size of y and preallocate memory
dim = size(y);
SingleInt = zeros(dim);
yFit = zeros(dim);
z = zeros(dim(2),1);
FWHMGauss = z; FWHMLorentz = z; x0 = z; IArea = z;

DGammaG = z; DGammaL = z; Dx0 = z; DeltaIArea = z;

for i=1:dim(2)
    fitresult = VoigtFit(x,y(:,i),'B0', 3348,'Gamma',2);
    % extracts parameters from fitresults
    yFit(:,i)=feval(fitresult,x);
    p=num2cell(abs(coeffvalues(fitresult)));
    [FWHMGauss(i),FWHMLorentz(i),a,x0(i)]=deal(p{:});

    % calculates first integral
    SingleInt(:,i)=a*Voigt(x,x0(i),FWHMGauss(i),FWHMLorentz(i),0);

    % calculates Area under first integral
    IArea(i)=a;

    % calculates uncertainties for fit parameters
    c=confint(fitresult);erros=abs(c(2,:)-c(1,:))/2;
    c=num2cell(erros);
    [DGammaG(i),DGammaL(i),DeltaIArea(i),Dx0(i)]=deal(c{:});
end

FitParam.x0 = x0;
FitParam.FWHMGauss = FWHMGauss;
FitParam.FWHMLorentz = FWHMLorentz;
FitParam.IArea = IArea;

FitError.Dx0 = Dx0;
FitError.DGammaG = DGammaG;
FitError.DGammaL = DGammaL;
FitError.DeltaIArea = DeltaIArea;

end