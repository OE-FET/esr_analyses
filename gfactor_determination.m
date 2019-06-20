function [g_sample, B_res]= gfactor_determination(x, y, pars, varargin)
%GFACTOR_DETERMINATION determines the g-factor of an ESR spectrum
%
%   [g_sample, B_res]= GFACTOR_DETERMINATION(x, y, pars)
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
%   SYNTAX:
%   [g_sample, B_res] = GFACTOR_DETERMINATION(x, y, pars)
%   [g_sample, B_res] = GFACTOR_DETERMINATION(x, y, pars, 'plot', 'y')
%
%   INPUT(S):
%   x        - Magnetic field in Gauss
%   y        - Derivative ESR spetrum
%   pars     - Structure containing measurement parameters
%   plot     - If plot == 'y', the spectrum is plotted with the resonance 
%              centre marked prominently.
%
%   OUTPUT(S):
%   g_sample - Sample g-factor
%   B_res    - Resonance center used to determine g-factor
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2019/05/06 12:58 $    $Revision: 1.1 $

plotting = get_kwarg(varargin, 'plot', 'n');

%% Resonance center from maximum of 1st integral as first guess
% integrate spectrum
[~, Int1] = double_int_num(x, y, 'baseline', 'n');
% find resonance peak at maximum
[~, II] = max(Int1,[],1);
B_res   = x(II);

%% g-factor from zero crossing, higher accuracy
% look at 0.4 Gauss interval around B_sample
dim = size(y);
pm = max(0.4, 5*pars.XSTEP);
Interval= [B_res-pm, B_res+pm];

for k=1:dim(2)
    Slice   = logical((Interval(k,1)<x) .* (x<Interval(k,2)));
    x_slice = x(Slice);
    y_slice = y(Slice,k);

    % find last value > 0 (B1)
    B1 = x_slice(sum(y_slice>0));
    % find first value < 0 (B2)
    B2 = x_slice(sum(y_slice>0)+1);

    % calculate intersect with x-axis
    I1 = y_slice(sum(y_slice>0));
    I2 = y_slice(sum(y_slice>0)+1);

    m = (I2-I1) / (B2-B1);
    c = I1 - m*B1;
    B_res(k) = -c/m;
end

% convert to g-factor
g_sample = b2g(B_res*1e-4, pars.MWFQ);

if plotting == 'y'
    for k = 1:dim(2)
        disp(['sample g = ', num2str(g_sample(k)), ', g_shift = ', ...
            num2str(round((g_sample(k)-gfree)*1e6)), ' ppm']);
    end
end

% plot the result
if plotting == 'y'
    figure(); hold on;
    [~, yoffsets] = stack_plot(x, y);
    xL = xlim; yL = ylim;
    stack_plot(B_res', zeros(length(B_res), 1)', 'yoffsets', yoffsets, 'style', 'ro');
    stack_plot(x, zeros(size(y)), 'yoffsets', yoffsets, 'style', 'k-');
    xlim(gca, xL); ylim(gca, yL);
    hold off;
end

end