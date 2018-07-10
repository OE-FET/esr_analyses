function [] = ReplaceTitleString(filePath)
%% [] = ReplaceQValue(Path2Data, Path2QVal)
% Replaces the dataset title string in the .DSC parameter file with the
% file title.
%
% Input arguments:
% Path2Data - .DTA/.DSC files with measurement data (saved by Xepr)
%

%% Read the .DSC parameter file

[folder, fileName] = fileparts(filePath);

Path2Par = [folder '/' fileName '.DSC'];
Path2Data = [folder '/' fileName '.DTA'];

[~, ~, Pars] = BrukerRead(Path2Data);

%% do nothing if Pars.TITL equals file name
if strcmp(Pars.TITL, sprintf('''%s''', fileName))
    return
end

%% replace Pars.TITL with file name
% read parameter file as text
A = regexp( fileread(Path2Par), '\n', 'split');
[~, d] = size(A);

% find TITL entry
for i = 1:d
    pos = strfind(A{i}, num2str(Pars.TITL));
    if pos > 0
        I = i;
    end
end
% replace current TITL with the file name
A{I} = sprintf('TITL	''%s''', fileName);

%% write changes to DSC file
fid = fopen(Path2Par, 'w');
fprintf(fid, '%s\n', A{:});
fclose(fid);

end