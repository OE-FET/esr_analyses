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
%           plotted.
%
%   OUTPUT(S):
%   s0          - In-phase signal from lock-in.
%   s90         - Out-of-phase-phase signal from lock-in.
%   phase_shift - Phase shift in rad used to calculate the above signals.
%

import esr_analyses.utils.*

phase_shift = get_kwarg(varargin, 'phase', false);
plot_result = get_kwarg(varargin, 'plot', false);

if ~phase_shift

    s90 = @(phi) sig_x * sin(phi) + sig_y * cos(phi);
    sig90amp = @(phi) sum(abs(s90(phi)) - mean(s90(phi)));

    phase_shift = fminbnd(sig90amp, -pi, pi);
end

s0 = sig_x * cos(phase_shift) + sig_y * sin(phase_shift);
s90 = sig_x * sin(phase_shift) + sig_y * cos(phase_shift);

if plot_result
    x = 1:length(sig_x);

    subplot(2, 1, 1)
    plot(x, sig_x, x, sig_y)
    legend(["Channel X", "Channel Y"])
    grid on; axis tight;

    subplot(2, 1, 2)
    plot(x, s0, x, s90)
    legend(["0 deg", "90 deg"])
    grid on; axis tight;
end

end
