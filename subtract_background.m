function [x, y, pars] = subtract_background(varargin)
%SUBTRACT_BACKGROUND Subtracts background signal from sample signal
% 	If desired, the result is written to a new Bruker ESR file. Experimental 
% 	conditions from DSC files are compared and a warning is issued when 
% 	differences are detected. ESR data is normalised for MW power, reciever 
% 	gain, number of scans, time constant, and modulation amplitude (to Hm = 1 G). 
% 	Before subtracting, the background signal is shifted to compensate for an
% 	offset in MWFQ.
%
% 	SYNTAX:
% 	subtract_background()              - prompts user for signal & background paths
% 	subtract_background(SPath, BPath)  - signal path, background path
%
% 	DEPENDENCIES:
% 	BrukerRead.m
% 	normalise_spectrum.m
% 	compare_pars.m
% 	stackplot.m
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2018/07/05 12:58 $    $Revision: 0.1 $

%%
global Path

% load files, prompt user if no file paths are given
switch nargin
    case 0
        [SName, SPath] = uigetfile([Path, '*.DTA'], 'Select signal data');
        File1 = [SPath, SName];
        [BName, BPath] = uigetfile([SPath, '*.DTA'], 'Select background data');
        File2 = [BPath, BName];
        Path = SPath;
        if isequal(BName, 0) %|| isequal(directory, 0)
            Path = [];
            fprintf('No background selected.');
            return
        end
    case 2
        File1 = varargin{1};
        File2 = varargin{2};
    otherwise
        error('Only 0 or 2 inputs are accepted.');
end

[xS, yS, pars] = BrukerRead(File1);
[xB, yB, ParsB] = BrukerRead(File2);

% normalise, if not yet done
[xS, yS, pars] = normalise_spectrum(xS, yS, pars);
[xB, yB, ParsB] = normalise_spectrum(xB, yB, ParsB);

% rescale background for Q-value, modulation amplitude and mw power
% WARNIG: All parameters affect both the signal amplitude and shape.
% Therefore, be caucious when subtracting a backround signal with 
% significantly different parameters.
Q_ratio = pars.QValue/ParsB.QValue;
B0MA_ratio = pars.B0MA/ParsB.B0MA;
Bmw_ratio = sqrt(pars.MWPW)/sqrt(ParsB.MWPW);

if Q_ratio > 1.2 || Q_ratio < 0.8
    disp(['Warning: Q-values differ by more than 20%.' newline ...
          'The Q-value may impact the signal shape and' newline ...
          'may prevent a proper background subtraction.']);
end

if Bmw_ratio > 1.2 || Bmw_ratio < 0.8
    disp(['Warning: MW fields differ by more than 20%.' newline ...
          'The MW field may impact the signal shape and' newline ...
          'may prevent a proper background subtraction.']);
end

if B0MA_ratio > 1.2 || B0MA_ratio < 0.8
    disp(['Warning: modulation amplitudes differ by more than 20%.' newline ...
          'The modulation amplitude may impact the signal shape and' newline ...
          'may prevent a proper background subtraction.']);
end

yB = yB * Q_ratio * B0MA_ratio * Bmw_ratio;

%% Compare experimental conditions
nDiff = compare_pars(pars, ParsB);

if nDiff > 0
    str = input('Do you want to continue ([y]/n)?', 's');
    if strcmpi(str, 'n') == 1
        error('Aborted.');
    end
end

%% Subtract spectra

%Compute and correct for H - offset
B_offset = (pars.MWFQ - ParsB.MWFQ) * planck/(gfree*bmagn) * 1E4;
disp(['B_offset = ', num2str(B_offset)]);
Bstep = xB(2) - xB(1);
offsetInterval = round(B_offset/Bstep);

if offsetInterval<0
    y = yS(1:end + offsetInterval, :) - yB(1 - offsetInterval:end, :);
    x = xS(1:end + offsetInterval);
elseif offsetInterval>0
    y = yS(1 + offsetInterval:end, :) - yB(1:end - offsetInterval, :);
    x = xS(1 + offsetInterval:end);
elseif offsetInterval == 0
    y = yS - yB;
    x = xS;
end

%% plot results
figure(1);
subplot(2, 1, 1);
hold on
yoffset = max(max(yS))*0.5;
% plot background
sp1 = stackplot(xB + B_offset, yB, 'yoffset', yoffset, 'style', 'r');
% plot signal
sp2 = stackplot(xS, yS, 'yoffset', yoffset, 'style', 'b'); 
hold off;
legend([sp1(1) sp2(1)],{'background', 'signal'})

% plot signal minus background
subplot(2, 1, 2)
stackplot(x, y, 'style', 'b');

end