function [xSub, ySub, ParsS] = SubtractBackground(varargin)
%SUBTRACTBACKGROUND Subtracts background signal from sample signal
% 	If desired, the result is written to a new Bruker ESR file. Experimental 
% 	conditions from DSC files are compared and a warning is issued when 
% 	differences are detected. ESR data is normalised for MW power, reciever 
% 	gain, number of scans, time constant, and modulation amplitude (to Hm = 1 G). 
% 	Before subtracting, the background signal is shifted to compensate for an
% 	offset in MWFQ.
%
% 	INPUT(S):
% 	SUBTRACTBACKGROUND()              - prompts user for signal & background paths
% 	SUBTRACTBACKGROUND(SPath, BPath)  - signal path, background path
%
% 	OUTPUT(S):
% 	x_sub                             - magnetic field [gauss]
% 	y_sub                             - subtracted and normalized signal
%
% 	DEPENDENCIES:
% 	BrukerRead.m
% 	NormaliseSpectrum.m
% 	comparePars.m
% 	StackPlot.m
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2018/07/05 12:58 $    $Revision: 0.1 $

%%
global Path

% load files, prompt user if no file paths are given
switch nargin
    case 0
        [SName, SPath] = uigetfile([Path, '*.DTA'],'Select signal data');
        File1 = [SPath, SName];
        [BName, BPath] = uigetfile([SPath, '*.DTA'],'Select background data');
        File2 = [BPath, BName];
        Path = SPath;
        if isequal(BName,0) %|| isequal(directory,0)
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

[xS, yS, ParsS] = BrukerRead(File1);
[xB, yB, ParsB] = BrukerRead(File2);

% normalise, if not yet done
[xS, yS, ParsS] = NormaliseSpectrum(xS, yS, ParsS);
[xB, yB, ParsB] = NormaliseSpectrum(xB, yB, ParsB);

% rescale for Q-value, modulation amplitude and mw power
% WARNIG: All parameters affect both the signal amplitude and shape.
% Therefore, be caucious when subtracting a backround signal with 
% significantly different parameters.
yB = yB * ParsS.QValue/ParsB.QValue * ParsS.B0MA/ParsB.B0MA * sqrt(ParsS.MWPW)/sqrt(ParsB.MWPW);


%% Compare experimental conditions
nDiff = comparePars(ParsS, ParsB);

if nDiff > 0
    str = input('Do you want to continue ([y]/n)?', 's');
    if strcmpi(str, 'n') == 1
        error('Aborted.');
    end
end

%% Subtract spectra

%Compute and correct for H - offset
B_offset = (ParsS.MWFQ - ParsB.MWFQ) * planck/(gfree*bmagn) * 1E4;
disp(['B_offset = ', num2str(B_offset)]);
Bstep = xB(2) - xB(1);
offsetInterval = round(B_offset/Bstep);

if offsetInterval<0
    ySub = yS(1:end + offsetInterval, :) - yB(1 - offsetInterval:end, :);
    xSub = xS(1:end + offsetInterval);
elseif offsetInterval>0
    ySub = yS(1 + offsetInterval:end, :) - yB(1:end - offsetInterval, :);
    xSub = xS(1 + offsetInterval:end);
elseif offsetInterval == 0
    ySub = yS - yB;
    xSub = xS;
end

%% plot results
hold off;
figure(1);

subplot(2, 1, 1);
yoffset = max(max(yS))*0.5;
% plot background
StackPlot(xB + B_offset, yB, 'yoffset', yoffset, 'style', 'r');
% plot signal
hold on; StackPlot(xS, yS, 'yoffset', yoffset, 'style', 'b'); hold off;
legend('background', 'signal')

% plot signal minus background
subplot(2, 1, 2)
StackPlot(xSub, ySub, 'style', 'b');

end