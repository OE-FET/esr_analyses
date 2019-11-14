function phase_cycle_gui(x, sig_x, sig_y)
%PHASE_CYCLE_GUI
%
%   Takes x- and y-channels from a lock-in amplifier and cycles them to
%   retrieve the in-phase and out-of-phase signals. Cycling is done
%   manually through a gui.
%
%   SYNTAX:
%   [s0, s90, phase_shift] = phase_cycle_gui(sig_x, sig_y)
%
%   INPUT(S):
%   sig_x - X-channel from lock-in.
%   sig_y - Y-channel from lock-in.
%

import esr_analyses.*
import esr_analyses.utils.*

s0 = @(phi) sig_x * cos(phi) + sig_y * sin(phi);
s90 = @(phi) sig_x * sin(phi) + sig_y * cos(phi);

sig0amp = @(phi) sum(abs(s0(phi) - mean(s0(phi))));
sig90amp = @(phi) sum(abs(s90(phi) - mean(s90(phi))));

phi_axis = -pi:2*pi/1000:pi;
phi_axis_deg = phi_axis*180/pi;
s0matrix = s0(phi_axis);
s90matrix = s90(phi_axis);

matrix = s0matrix; matrix(:, :, 2) = s90matrix;

subplot(2, 1, 1)
plot(phi_axis_deg, sig0amp(phi_axis), phi_axis_deg, sig90amp(phi_axis))
legend(["0 deg signal", "90 deg signal"])
axis tight
xlabel('Phase shift [deg]')
xticks(0:45:360)
ylabel('Signal strength')
subplot(2, 1, 2)

slider_plot(x, matrix, 'slider_axis', phi_axis_deg, 'slider_unit', ' deg');

end
