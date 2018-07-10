function [argout, SortedDATA, fitresult] = PowerSatAnalysesSliceFIT(varargin)
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
% Dependencies:
% BrukerRead.m
% NormaliseSpectrum.m
% SubtractBackground.m
% DoubleIntFIT.m
% gfactor_determination.m
% MWmean.m

%% Read data 
if nargin == 1
    Path = varargin{1};
    str = input('Would you like to subtract a background signal? y/[n]','s');
    if strcmp(str,'y')==1
        [x,y,Pars]=SubtractBackground(Path);
    else
        [x,y,Pars]=NormaliseSpectrum(Path);
    end
elseif nargin == 0
    str = input('Would you like to subtract a background signal? y/[n]','s');
    if strcmp(str,'y')==1
    [x,y,Pars]=SubtractBackground;
    else
    [x,y,Pars]=NormaliseSpectrum;
    end
elseif nargin ==3
    x = varargin{1}; y = varargin{2}; Pars = varargin{3};
    [x,y,Pars]=NormaliseSpectrum(x,y,Pars);
end

% Get Q value powers
try
    QValue = Pars.QValue;
catch
    QValue = input('Please give cavity QValue:');
end

% Get MW powers
if strcmp(Pars.YTYP,'IGD')==1
    Pmw = Pars.z_axis * 1E-3; % MW Power in W
else
    error('The specified file is not a 2D data file.');
end

%% calculate integrated intensities

% normalise for measurement conditions
[x,y] = NormaliseSpectrum(x,y,Pars);

% double integration of spectra
fitresult = PseudoVoigtFit(x,y);
for i = 1:size(y,2)
    Ipp(i) = fitresult{i}.a; % peak to peak amplitude
    s(i) = fitresult{i}.s; % gaussness
    Hpp(i) = fitresult{i}.w; % peak to peak line width
    H0(i) = fitresult{i}.x0; % resonance center in Gauss
    DoubleInt(i) = (1-s(i)) * ( Ipp(i)*pi*Hpp(i)^2/sqrt(3) ) + s(i) * ( Ipp(i)*exp(1/2)*(1/4)*sqrt(pi/2)*Hpp(i)^2 );

    error_matrix = confint(fitresult{i},0.95);
    errors = (error_matrix(2,:)-error_matrix(1,:))/2;
    
    dIpp(i) = errors(1);
    ds(i) = errors(2);
    dHpp(i) = errors(3);
    dH0(i) = errors(4);
    dDoubleInt(i) = dIpp(i)*pi*Hpp(i)^2/sqrt(3) + dHpp(i)*2*Ipp(i)*pi*Hpp(i)/sqrt(3); % approximate
end

% determine g-value
g = b2g(H0*1e-4, Pars.MWFQ);
gDelta = g/(H0.*1e-4).* dH0;

% convert MWPW to magnetic field strength in Tesla
% use Qref = 7500 and c = 2.0 for Nagoya files!!!
Qref = 8355;
c = 2.2;
cCryo = 2.2*1.3; % calibrated with gMarker

f_mean = MWmean(Pars);

Bmw = f_mean * cCryo * sqrt(Pmw) * sqrt(QValue/Qref) * 1E-4; % in Tesla

% collect values in output matrix for easy copying to Origin
output = [Pmw Bmw DoubleInt' dDoubleInt' Ipp' dIpp' Hpp' dHpp' g' gDelta'];

xNorm = x;
yNorm = y;


%% combine data from several files

% DATA = cell2mat(output);
DATA = output;

% sort data with ascending MW Power

[~,I]=sort(DATA(:,1));
SortedDATA=DATA(I,:);

argout.MircowavePower = SortedDATA(:,1);
argout.MicrowaveB = SortedDATA(:,2);
argout.DoubleIntegrated = SortedDATA(:,3);
argout.PeakIntensity = SortedDATA(:,5);
argout.p2pLineWidth = SortedDATA(:,7);
argout.gValue = SortedDATA(:,9);
argout.xNorm = xNorm;
argout.yNorm = yNorm;

%% Fit and plot results
gValue = mean(argout.gValue);
try
fitresult = saturationCurveFit(argout.MicrowaveB(1:end),...
    abs(argout.DoubleIntegrated(1:end)), gValue);
catch
    disp('Fit could not be performed.');
end