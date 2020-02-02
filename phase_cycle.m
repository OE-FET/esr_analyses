function [s0, s90, phase_shift] = phase_cycle(sig_x, sig_y, varargin)
%PHASE_CYCLE
%
%   Takes x- and y-channels from a lock-in amplifier and cycles them by the
%   given phase to retrieve the in-phase and out-of-phase signals. If no
%   phase is given, the phase which minimizes the out-of-phase signal is
%   determined automatically.
%
%   SYNTAX:
%   [s0, s90, phase_shift] = phase_cycle(sig_x, sig_y)
%   [s0, s90, phase_shift] = phase_cycle(sig_x, sig_y, 'phase', phase)
%   [s0, s90, phase_shift] = phase_cycle(sig_x, sig_y, 'plot', true)
%
%   INPUT(S):
%   sig_x - X-channel from lock-in.
%   sig_y - Y-channel from lock-in.
%   phase - Phase offset in rad to use. If not given, the appropriate
%           phase will be determined by minimizing the out-of-phase signal.
%   plot  - If true, the original and phase-corrected signals will be
%           plotted. Defaults to true.
%
%   OUTPUT(S):
%   s0          - In-phase signal from lock-in.
%   s90         - Out-of-phase-phase signal from lock-in.
%   phase_shift - Phase shift in rad used to calculate the above signals.
%

import esr_analyses.*
import esr_analyses.utils.*

phase_shift = get_kwarg(varargin, 'phase', false);
plot_result = get_kwarg(varargin, 'plot', true);

if ~phase_shift
    s90 = @(phi) sig_x * sin(phi) + sig_y * cos(phi);
    sig90amp = @(phi) sum(s90(phi).^2, 'all');

    phase_shift = fminbnd(sig90amp, -pi, pi) + pi;
end

s0 = sig_x * cos(phase_shift) - sig_y * sin(phase_shift);
s90 = sig_x * sin(phase_shift) + sig_y * cos(phase_shift);

if plot_result
    
    fig_name = 'Phase cycling';
    
    fh = findobj( 'Type', 'Figure', 'Name', fig_name);
    if isempty(fh)
        figure('Name', fig_name);
    end
    
    x = (1:length(sig_x))';
    
    yoffsetsX = max(sig_x);
    yoffsetsY = max(sig_y);
    yoffsets = max([yoffsetsX; yoffsetsY])*2;
    
    hold off;

    subplot(1, 2, 1)
    h1 = stackplot(x, sig_x, 'yoffsets', yoffsets, 'style', 'b');
    hold on;
    h2 = stackplot(x, sig_y, 'yoffsets', yoffsets, 'style', 'r');
    legend([h1(1), h2(1)], ["Channel X", "Channel Y"])
    grid on; axis tight;

    subplot(1, 2, 2)
    h1 = stackplot(x, s0, 'yoffsets', yoffsets, 'style', 'b');
    hold on;
    h2 = stackplot(x, s90, 'yoffsets', yoffsets, 'style', 'r');
    legend([h1(1), h2(1)], ["0 deg", "90 deg"])
    grid on; axis tight;
    
    sgtitle('Phase cycling')
    
    hold off;
    
    rad = input('Apply correction? [0 deg] ');
    
    if ~isempty(rad)
        rad = rad*pi/180;
        [s0, s90, phase_shift] = phase_cycle(sig_x, sig_y, 'phase', phase_shift+rad, 'plot', true);
    end
end

phase_shift = mod(phase_shift, 2*pi);

end
