function [x, y, parsS] = subtract_background(varargin)
%SUBTRACT_BACKGROUND Subtracts background signal from sample signal.
%
% 	If desired, the result is written to a new Bruker ESR file. Experimental
% 	conditions from DSC files are compared and a warning is issued when
% 	differences are detected. ESR data are normalised for MW power, reciever
% 	gain, number of scans, time constant, and modulation amplitude (to Hm = 1 G).
% 	Before subtracting, the background signal is shifted to compensate for an
% 	offset in MWFQ.
%
% 	INPUT(S):
%   SUBTRACT_BACKGROUND()           - opens gui to select signal and bg paths
%   ...('signal_data')              - given signal path & opens a gui to
%                                     select bg path
%   ...('signal_data','bg_data')    - given signal & bg paths
%
%   OUTPUT(S):
% 	x       - magnetic field axis
%   y       - normalised & subtracted signal intensity
%   pars    - experimental parameters
%
%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2019/05/06 12:58 $    $Revision: 1.1 $

import esr_analyses.*
import esr_analyses.utils.*

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
    case 1
        File1 = varargin{1};
        [SPath,~,~] = fileparts(File1);
        [BName, BPath] = uigetfile([SPath, '*.DTA'], 'Select background data');
        File2 = [BPath, BName];
    case 2
        File1 = varargin{1};
        File2 = varargin{2};
    otherwise
        error('Only 0-2 inputs are accepted.');
end

[xS, yS, parsS] = BrukerRead(File1);
[xB, yB, ParsB] = BrukerRead(File2);

% normalise, if not yet done
[xS, yS, parsS] = normalise_spectrum(xS, yS, parsS);
[xB, yB, ParsB] = normalise_spectrum(xB, yB, ParsB);

% rescale background for Q-value, modulation amplitude and mw power
% WARNING: All parameters affect both the signal amplitude and shape.
% Therefore, be cautious when subtracting a backround signal with
% significantly different parameters.
Q_ratio = parsS.QValue/ParsB.QValue;
B0MA_ratio = parsS.B0MA/ParsB.B0MA;
Bmw_ratio = sqrt(parsS.MWPW)/sqrt(ParsB.MWPW);
ratiosary = [Q_ratio, B0MA_ratio, Bmw_ratio];
ratiostxt = {'Q-values','Modulation amplitudes','Microwave fields'};
for ii = 1:length(ratiosary)
    if ratiosary(ii) > 1.2 || ratiosary(ii) < 0.8
        warning([ratiostxt{ii} ' differ by more than 20%. '...
                 'This may impact the signal shape and '...
                 'prevent a proper background subtraction.']);
    end
end

yB = yB * Q_ratio * B0MA_ratio * Bmw_ratio;

%% Compare experimental conditions
nDiff = compare_pars(parsS, ParsB);

if nDiff > 0
    str = input('Do you want to continue ([y]/n)?', 's');
    if strcmpi(str, 'n')
        error('Aborted.');
    end
end

%% Subtract spectra

% Compute and correct for the expected resonance centre offset due to
% differing microwave fields (from ge*bmagn*B = hbar*f)
B_offset = (parsS.MWFQ - ParsB.MWFQ) * planck/(gfree*bmagn) * 1E4;
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
figure('Name','Background subtraction');
subplot(2, 1, 1);
hold on
% 27/07/19: For PowerSat data: It appears as if you meant to use the below
% line to send the offsets to stackplot, but you used the keyword "yoffset"
% rather than "yoffsets" so the offset wasn't interpreted in stackplot and
% instead defaulted to stackplot's calculation of yoffsets = [0
% max(ydiff)*1.3]. This caused different offsets for the bg and signal
% data (and the slices didn't align). Moreoever, max(max(yS)) was returning
% one value, meaning that even if you sent it properly to stackplot, it
% was sending one value rather than an array. Every slice would be moved up
% by the same amount (and thus they would not be offset). Hopefully this
% change is what you meant to do: every slice is moved up by a different
% amount and the bg and signal slices overlap.
yoffset = max(yS)*0.5;
% plot background
sp1 = stackplot(xB + B_offset, yB, 'yoffsets', yoffset, 'style', 'r');
% plot signal
sp2 = stackplot(xS, yS, 'yoffsets', yoffset, 'style', 'b');
hold off;
legend([sp1(1) sp2(1)],{'background', 'signal'})

% plot signal minus background
subplot(2, 1, 2)
stackplot(x, y, 'style', 'b');

end