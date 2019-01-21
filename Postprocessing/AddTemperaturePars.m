function [] = AddTemperaturePars(filePath)
%% [] = ReplaceQValue(Path2Data, Path2QVal)
% Reads the measurement temperature from the file name as '000K' and adds it
% as Cryostat information layer to the end of the .DSC file.
%
% Input arguments:
% Path2Data - .DTA/.DSC files with measurement data (saved by Xepr)
%

%% Get measurement temperature from file name or .DSC entry

[folder, fileName] = fileparts(filePath);

Path2Par = [folder '/' fileName '.DSC'];
Path2Data = [folder '/' fileName '.DTA'];

[~, ~, Pars] = BrukerRead(Path2Data);

if isfield(Pars, 'Temperature')
    return
else
    [startIndex, endIndex] = regexp(Path2Data, '[0-9]+K');
    T = str2double( Path2Data(startIndex:endIndex-1) );
    if isnan(T)
        T = 298;
    end
    Pars.Temperature = [num2str(T), ' K'];
end

% read parameter file
A = regexp( fileread(Path2Par), '\n', 'split');

%% Append temperature data if not present

Iend = length(A);

A{Iend-2} = '.DVC     tempCtrl, 1.0'                                        ;
A{Iend-1} = ''                                                              ;
A{Iend  } = 'AcqWaitTime        120.0 s'                                    ;
A{Iend+1} = sprintf('Temperature        %.1f K', T)                         ;
A{Iend+2} = 'Tolerance          0.1 K'                                      ;
A{Iend+3} = ''                                                              ;
A{Iend+4} = '*'                                                             ;
A{Iend+5} = '************************************************************'  ;


%% write changes to DSC file
fid = fopen(Path2Par, 'w');
fprintf(fid, '%s\n', A{:});
fclose(fid);

end