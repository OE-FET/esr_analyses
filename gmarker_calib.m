function [xnew, y, pars] = gmarker_calib(x, y, pars, varargin)
%GMARKER_CALIB Calibrates the magnetic field axis with respect to a marker.
%
%   If the optional input 'gMarker' is not given, a standard Bruker marker
%   with a g-value of 1.979843 is assumed.
%
% 	SYNTAX:
% 	[xnew, y, pars] = GMARKER_CALIB(x, y, pars)
% 	[xnew, y, pars] = GMARKER_CALIB(x, y, pars, 'gMarker', marker_g_value)
%
% 	OUTPUT(S):
% 	xnew		- magnetic field axis normalized by g-marker
%	y			- ESR spectrum
%	pars		- aquisition parameter structure
%	                                  

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2019/05/06 12:58 $    $Revision: 1.1 $

%% Input analyses
marker_g = get_kwarg(varargin, 'gMarker', 1.979843);

%% Find resonance peak of marker
B_target = 10E3 * g2b(marker_g,pars.MWFQ);

% find 10 Gauss interval around suspected marker location
Interval = [B_target-5, B_target+5];
Slice = logical((Interval(1)<x) .* (x<Interval(2)));
x_slice = x(Slice);
y_slice = y(Slice);

if isempty(x_slice)
    disp('Marker location outside of spectrum. Proceeding without marker calibration.');
    xnew=x;
    return;
end

% heuristic check if there is a marker ESR signal
y_slice_int = cumtrapz(x_slice, y_slice);
y_slice_int_norm = y_slice_int/max(y_slice_int);
[~, ~, ~, p] = findpeaks(y_slice_int_norm);
if p < 0.3
    disp('No ESR signal found at marker location.');
    plot(x, y, x_slice, y_slice);
    xnew=x;
    return;
end

% determine apparent g-factor of marker
[~, B_marker] = gfactor_determination(x_slice, y_slice, pars);

%% perform shift of x-axis
% calculate x-axis offset
B_offset = B_target - B_marker;
disp(['B_offset = ', num2str(B_offset), 'G']);

% save new x-axis, confirm
xnew = x + B_offset;
disp('Offset corrected.');

% plot spectrum with marker location for visual confirmation
hold off; plot(xnew, y, 'b');
hold on; plot(x, y, 'b--');
plot(B_target, 0, 'o'); plot(B_marker ,0, 'o');
plot(xnew, zeros(size(xnew)), 'k','LineWidth', 1);
axis tight;
legend('Corrected spectrum', 'Original spectrum', ...
    'Corrected marker resonance', 'Original marker resonance');
hold off;

end