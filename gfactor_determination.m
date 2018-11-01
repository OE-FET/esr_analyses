function [g_sample, B_res]= gfactor_determination(x, y, Pars, plotting)
%GFACTOR_DETERMINATION determines the g-factor of an ESR spectrum
%   [g_sample, B_res]= GFACTOR_DETERMINATION(x, y, Pars, plotting)
%   determines the g-factor of derivative ESR spectra given by (x,y) with
%   the microwave frequency in the paramater structure PARS. G_SAMPLE is
%   the g-factor, B_RES is the resonance field.
%
%   The resonance center is detrmined intially by finding the maximum of the
%   integrated spectrum and is refinded by finding the closest zero-crossing
%   of the derivative spectrum. The location of B_RES is interpolated
%   from the two closest data-points with x > 0 and x < 0. The resonance
%   center field is converted to a g-factor with the function B2B.
%
%   This function also works with hyperfine-split derivative spectra that
%   have multiple zero-crossings.
%
%   INPUT:
%   X        - Magnetic field in Gauss.
%   Y        - Derivative ESR spetrum.
%   PARS     - Structure containing measurement parameters. It must contain 
%              the microwave frequency Pars.MWFQ in Hz.
%   VARARGIN - 'y' or 'n' determines if the results are plotted.
%              Default is 'n'.
%
%   OUTPUT:
%   G_FACTOR - Sample g-factor.
%   B_RES    - Resonance center used to determine g-factor.
%
%
%   DEPENDENCIES:
%   StackPlot.m
%   DoubleIntNUM.m
%   b2g.m
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2018/07/05 12:58 $    $Revision: 0.1 $

if nargin < 4
    plotting = 'y';
end

%% Resonance center from maximum of 1st integral as first guess
% integrate spectrum
[~, Int1] = DoubleIntNUM(x, y, 'n');
% find resonance peak at maximum
[~, II] = max(Int1,[],1);
B_res=x(II);

%% g-factor from zero crossing, higher accuracy
% look at 0.4 Gauss interval around B_sample
dim = size(y);
pm = 0.4;
Interval= [B_res-pm, B_res+pm];

for k=1:dim(2)
    Slice = logical((Interval(k,1)<x) .* (x<Interval(k,2)));
    x_slice = x(Slice);
    y_slice = y(Slice,k);

    % find last value > 0 (B1)
    B1 = x_slice(sum(y_slice>0));
    % find first value < 0 (B2)
    B2 = x_slice(sum(y_slice>0)+1);

    % determine slope m for line between B1 and B2
    I1=y_slice(sum(y_slice>0));
    I2=y_slice(sum(y_slice>0)+1);
    m = (I2-I1)/(B2-B1); c=I1-m*B1;
    % calculate intersect with x-axis
    B_res(k)=-c/m;
end
% convert to g-factor
g_sample=b2g(B_res*1e-4,Pars.MWFQ);

if plotting=='y'
    for k=1:dim(2)
        disp(['sample g = ',num2str(g_sample(k)),', g_shift = ',num2str(round((g_sample(k)-gfree)*1e6)),' ppm']);
    end
end

% plot the result
if plotting=='y'
    yoffset = 0.5*max(max(y));
    hold off; StackPlot(x,y,'yoffset',yoffset);hold on;
    xL=xlim;yL=ylim;
    [~,h1]=StackPlot(B_res',zeros(length(B_res),1)','yoffset',yoffset);
    set(h1,'Marker','o','Color','r');
    [~,h2]=StackPlot(x,zeros(size(y)),'yoffset',yoffset);
    set(h2,'Color','k','LineWidth',1);
    xlim(gca,xL);ylim(gca,yL);
    hold off;
end

end