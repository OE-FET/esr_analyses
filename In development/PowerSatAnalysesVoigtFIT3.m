function [argout] = PowerSatAnalysesVoigtFIT3(varargin)
%% Power saturation curves, analytical integration form fit
% Calculates power stauration curves of integrated intensity and maximum
% intensity. This programm fits the ESR signal with a Voigt curve derivative
% and then performs analytical integration.
% Advantage: Long tails of the resonance peak are not negletcted.
% Disadvantage: Resonance curve has to be symmetrical and of a Voigt-shape.
%
% Power saturation curves are given in SortedDATA with the colmuns:
% [Pmw MwB II I_max Hpp g]
%
% All spectra are collected in the structure 'argout.ERSIntensity'.
%
% DEPENDENCIES:
% BrukerRead.m
% NormaliseSpectrum.m
% SubtractBackground.m
% ESRVoigtSimulation.m
% gfactor_determination.m
% MWmean.m
% getSamplePosition.m
% Natural Constants

%% Read data 
if nargin == 1
    Path = varargin{1};
    str = input('Would you like to subtract a background signal? y/[n]','s');
    if strcmp(str,'y')==1
        [x, y, Pars] = SubtractBackground(Path);
    else
        [x, y, Pars] = NormaliseSpectrum(Path);
    end
elseif nargin == 0
    str = input('Would you like to subtract a background signal? y/[n]','s');
    if strcmp(str,'y')==1
    [x, y, Pars] = SubtractBackground;
    else
    [x, y, Pars] = NormaliseSpectrum;
    end
elseif nargin ==3
    x = varargin{1}; y = varargin{2}; Pars = varargin{3};
    [x, y, Pars] = NormaliseSpectrum(x, y, Pars);
end

%% Get Q value
try
    QValue = Pars.QValue;
catch
    QValue = input('Please give cavity QValue:');
end

%% Get MW powers
if strcmp(Pars.YTYP, 'IGD')==1
    Pmw = Pars.z_axis/1000; % MW Power in W
else
    error('The specified file is not a 2D data file.');
end

% ask for height and length of sample if not in Pars

Pars = getSamplePosition(Pars);

% convert MWPW to magnetic field strength in Tesla
% use Qref = 7500 and c = 2.0 for Nagoya files
Qref = 8355;
c = 2.2;

position_correction = MWmean(Pars);

Bmw = position_correction * c * sqrt(Pmw) * sqrt(QValue/Qref) * 1E-4; % in Tesla

%% Get staring points for fit

% normalise for measurement conditions
[x, y] = NormaliseSpectrum(x, y, Pars);

% get initial parameters for fit
mid = round(length(Bmw)/2);
fit1 = PseudoVoigtFit(x, y(:,mid));

Ipp = abs(fit1.a); % peak to peak amplitude in Gauss
Hpp = fit1.w; % peak to peak line width in Gauss
T2 = 2/sqrt(3) * 1/(gmratio * Hpp*1E-4); % starting point for fit, sec
T1 = 2*T2; % starting point for fit, sec
B0 = fit1.x0; % resonance center in Gauss

[X, Y] = meshgrid(x, Bmw); Z = y; % grid data for fitting algorithm


%% Perform Voigt fit
A = [0.0103 0.0103 0.0150]*2/3;
B0 = [3343.31293941321 3343.31293941321 3343.36693724911];
T1 = [6.17196041282668e-06 6.17196041282668e-06	0.000133918005979904];
T2 = [1.34901439609167e-07 1.34901439609167e-07	1.98816404899825e-07];
Brms = [0.2 0.202292914879847	0.296353090180823];

var0 = [A(1) B0(1) T1(1) T2(1) Brms(1) A(2) B0(2) T1(2) T2(2) Brms(2) A(3) B0(3) T1(3) T2(3) Brms(3)]; % vector with starting points [A, T2, Brms]
var0mod = var0 + (rand(size(var0))-0.5).*var0*1.5e-4;

% function to minimize: sum of squared errors, two peaks
fitfunc = @(var) abs(var(1))*ESRVoigtSimulation(X , var(2),  abs(var(3)),  abs(var(4)), Y, abs(var(5)), 1, Pars.B0MA*1e4) + abs(var(6))*ESRVoigtSimulation(X , var(7), abs(var(8)), abs(var(9)), Y, abs(var(10)), 1, Pars.B0MA*1e4) + abs(var(11))*ESRVoigtSimulation(X , var(12), abs(var(13)), abs(var(14)), Y, abs(var(15)), 1, Pars.B0MA*1e4);
sumofsquares = @(var) sum(sum( abs(fitfunc(var) - Z).^2  ));

% Fit model to data with fminsearch (Nelder Mead algorithm, much better
% convergance than Levenberg Marquard or trust Region)
opt = optimset('TolFun',1e-8,'TolX',1e-8,'PlotFcns',@optimplotfval, 'MaxFunEvals', 1e7, 'MaxIter', 1e7);
[ft_rslt, sumofsquares_error] = fminsearch(sumofsquares, var0, opt);

A = [abs(ft_rslt(1)), abs(ft_rslt(6)), abs(ft_rslt(11))];
B0 = [abs(ft_rslt(2)), abs(ft_rslt(7)), abs(ft_rslt(12))];
T1 = [abs(ft_rslt(3)), abs(ft_rslt(8)), abs(ft_rslt(13))];
T2 = [abs(ft_rslt(4)), abs(ft_rslt(9)), abs(ft_rslt(14))];
Brms = [abs(ft_rslt(5)), abs(ft_rslt(10)), abs(ft_rslt(15))];

% get best fit curve
zFit = fitfunc(ft_rslt);
for i=1:length(var0)/5
    name = sprintf('comp%i',i);
    zFitComponent.(name) = A(i)*ESRVoigtSimulation(X , B0(i), T1(i), T2(i), Y, Brms(i), 1, Pars.B0MA*1e4);
end

%% Estimate fitting errors
dof = length(x) - length(var0); % degrees of freedom in fitting problem
sdr = sqrt(sumofsquares_error/dof); % standard deviation of residuals
J = jacobianest(fitfunc, ft_rslt); % jacobian matrix
Sigma = sdr^2*inv(J'*J); % covariance matrix
se = sqrt(diag(Sigma))'; % parameter standrad errors

%% plot results
Xt = X'; Yt = Y';

f1 = figure( 'Name', '3D Voigt Fit' );
figure(f1);
hold off;
scatter3(Xt(:), Yt(:), Z(:),'.k');
hold on;
surf(X, Y, zFit','FaceAlpha', 0.5, 'EdgeColor','none');
legend('Data', 'Fit','Location','northeast')
xlabel('Magnetic Field [G]')
ylabel('Microwave Field Field [T]')
zlabel('ESR signal [a.u.]')
axis tight;

f2 = figure( 'Name', '3D Voigt Fit, x-section' );
figure(f2);
offset = max(max(y))*0.8;
hold off;
StackPlot(x,y,'yoffset',offset,'style','.k');
hold on;
StackPlot(x,zFit,'yoffset',offset,'style','-r');
StackPlot(x,zFitComponent.comp2,'yoffset',offset,'style','-b');
StackPlot(x,zFitComponent.comp3,'yoffset',offset,'style','-g');
legend('Data');

%% Susceptebility Calculation
for i=1:length(var0)/5
    FWHMLorentz = 2/(gmratio*T2(i)) * sqrt(1 + gmratio^2*Bmw.^2*T1(i)*T2(i)); % in Tesla
    s = 1246.572123192064; % scaling factor for pseudo modulation
    FittedArea(:,i) = s * Pars.B0MA * 1e4 * Bmw./FWHMLorentz .* A(i);
    Pars.GFactor = b2g(B0(i)*1e-4, Pars.MWFQ);

    [ChiAll, dChiAll] = SusceptebilityCalc(FittedArea(:,i), Pars);
    [NSpinAll, dNSpinAll] = SpinCounting(FittedArea(:,i), Pars);
    
    Chi(i) = ChiAll(1);
    dChi(i) = dChiAll(1);
    NSpin(i) = NSpinAll(1);
    dNSpin(i) = dNSpinAll(1);
end

%% Output

% create output structure
argout.x = x;
argout.y = y;
argout.yFit = zFit;
argout.components = zFitComponent;
argout.Pars = Pars;

argout.A = A;
argout.T1 = T1;
argout.T2 = T2;
argout.Brms = Brms;
argout.B0 = B0;

argout.dA = [se(1), se(5)];
argout.dT1 = [se(2), se(6)];
argout.dT2 = [se(3), se(7)];
argout.dBrms = [se(4), se(8)];
argout.dB0 = [se(9), se(10)];

argout.gfactor = Pars.GFactor;
argout.Bmw = Bmw;

argout.Chi = Chi;
argout.dChi = dChi;
argout.NSpin = NSpin;
argout.dNSpin = dNSpin;

end
