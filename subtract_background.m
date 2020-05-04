function dset = subtract_background(varargin)
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
%   SUBTRACT_BACKGROUND()          - opens GUI to select signal and
%                                    background data files
%   ...('signal_data')             - given path to signal data & opens a 
%                                    GUI to select background data file
%   ...('signal_data','bg_data')   - given signal & bg data files
%   ...(dsetSig, dsetBg)           - given signal and bg data directly
%
%   OUTPUT(S):
% 	dset - new dataset
%
%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$

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
        dsetSig = BrukerRead(File1);
        dsetBG = BrukerRead(File2);
    case 1
        File1 = varargin{1};
        [SPath,~,~] = fileparts(File1);
        [BName, BPath] = uigetfile(fullfile(SPath, '*.DTA'), 'Select background data');
        File2 = [BPath, BName];
        dsetSig = BrukerRead(File1);
        dsetBG = BrukerRead(File2);
    case 2
        if istable(varargin{1}) && istable(varargin{2})
            dsetSig = varargin{1};
            dsetBG = varargin{2};
        else
            File1 = varargin{1};
            File2 = varargin{2};
            dsetSig = BrukerRead(File1);
            dsetBG = BrukerRead(File2);
        end
    otherwise
        error('Please give either data sets or file paths as input.');
end

dsetSig = normalise_spectrum(dsetSig);
dsetBG = normalise_spectrum(dsetBG);

parsS = dsetSig.Properties.UserData;
parsB = dsetBG.Properties.UserData;

% rescale background for Q-value, modulation amplitude and mw power
% WARNING: All parameters affect both the signal amplitude and shape.
% Therefore, be cautious when subtracting a backround signal with
% significantly different parameters.
Q_ratio = parsS.QValue/parsB.QValue;
B0MA_ratio = parsS.B0MA/parsB.B0MA;
Bmw_ratio = sqrt(parsS.MWPW)/sqrt(parsB.MWPW);
ratiosary = [Q_ratio, B0MA_ratio, Bmw_ratio];
ratiostxt = {'Q-values','Modulation amplitudes','Microwave fields'};
for ii = 1:length(ratiosary)
    if ratiosary(ii) > 1.2 || ratiosary(ii) < 0.8
        warning([ratiostxt{ii} ' differ by more than 20%. '...
                 'This may impact the signal shape and '...
                 'prevent a proper background subtraction.']);
    end
end

dsetBG{:,2:end} = dsetBG{:,2:end} * Q_ratio * B0MA_ratio * Bmw_ratio;

%% Compare experimental conditions
compare_pars(parsS, parsB);

% if nDiff > 0
%     str = input('Do you want to continue ([y]/n)?', 's');
%     if strcmpi(str, 'n')
%         error('Aborted.');
%     end
% end

%% Subtract spectra

dset = dsetSig;

% Compute and correct for the expected resonance centre offset due to
% differing microwave fields (from ge*bmagn*B = hbar*f)
B_offset = (parsS.MWFQ - parsB.MWFQ) * planck/(gfree*bmagn) * 1E4;
disp(['B_offset = ', num2str(B_offset)]);
Bstep = dsetSig{1,1} - dsetSig{2,1};
offsetInterval = round(B_offset/Bstep);

if offsetInterval < 0
    y = dsetSig{:,2:end}(1:end + offsetInterval, :) - dsetBG{:,2:end}(1 - offsetInterval:end, :);
    x = dsetSig{:,1}(1:end + offsetInterval, :);
elseif offsetInterval > 0
    y = dsetSig{:,2:end}(1 + offsetInterval:end, :) - dsetBG{:,2:end}(1:end - offsetInterval, :);
    x = dsetSig{:,1}(1 + offsetInterval:end);
elseif offsetInterval == 0
    y = dsetSig{:,2:end} - dsetBG{:,2:end};
    x = dsetSig{:,1};
end

n = length(x);

dset{1:n,:} = [x, y];
dset(n+1:end, :) = [];

%% plot results

for k=1:width(dset)-1
    figure('Name','Background subtraction');
    
    subplot(2, k, k)

    yoffset = max(dsetSig{:,k+1});
    % plot background
    sp1 = stackplot_xepr(dsetBG(:,[1, k+1]), 'yoffsets', yoffset, 'style', 'r');
    hold on;
    % plot signal
    sp2 = stackplot_xepr(gca, dsetSig(:,[1, k+1]), 'yoffsets', yoffset, 'style', 'b');
    hold off;
    legend([sp1(1,1) sp2(1,1)], {'background', 'signal'})
end

end