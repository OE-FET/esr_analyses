function [x_scaled, y_scaled, pars] = to_current_scale(x, y, pars, varargin)
% Rescales the y-axis from a.u. in Xepr to change in EDMR current [Ampere /
% Gauss]. Input signal should already be normalised.

import esr_analyses.*
import esr_analyses.utils.*

amps2volts = get_kwarg(varargin, 'amps2volts', 1e6);  % from transimpedance amplifier

y_volts = y * 1/1000;  % 1000 a.u. corresponds to 2 Vpp

y_amps = y_volts ./ amps2volts;
y_amps_per_gauss = y_amps ./ (pars.B0MA*1e4);

pars.IRNAM = 'dI/dB';
pars.IRUNI = 'A/G';

x_scaled = x;
y_scaled = y_amps_per_gauss;

end