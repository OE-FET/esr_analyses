function [xnew,y,Pars] = MarkerCalib(varargin)
%MARKERCALIB Calibrates the magnetic field axis with respect to a marker signal.
% If the optional input 'gMarker' is not given, a standard Bruker marker with a
% g-value of 1.979843 is assumed.
%
% 	INPUT(S):
% 	MARKERCALIB()                                  - prompts user for Xepr file
% 	MARKERCALIB('/path/to/file')                   - path to Xepr file
% 	ARKERCALIB(x,y,Pars)                           - B field, signal, params.
% 	MARKERCALIB(x,y,Pars,'gMarker',marker_g_value) - ..., & uses given g-value
%                                                  as marker g-value
%
% 	OUTPUT(S):
% 	xnew		- magnetic field axis normalized by g-marker
%	y			- ESR spectrum
%	Pars		- aquisition parameter structure
%	                                  
% 
% 	DEPENDENCIES:
% 	BrukerRead.m
% 	g2b.m
% 	DoubleIntNUM.mg
% 	factor_determination.m
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2018/07/05 12:58 $    $Revision: 0.1 $

%% Input analyses
marker_g = 1.979843; nInput = nargin;
% If gmarker is specified then set value.
intInd = find(strcmp('gMarker', varargin));
if isempty(intInd) == 0
    marker_g = varargin{intInd+1};
    nInput=nargin-2;
end

switch nInput
    case 0
        [x,y,Pars] = BrukerRead;
    case 1
        Path=varargin{1};
        [x,y,Pars] = BrukerRead(Path);
    case 3
        x = varargin{1};y = varargin{2};Pars = varargin{3};
    otherwise
        error(['Invalid arguments. Please give either path to filename or '...
               '(x,y,Pars) as input. ''gMarker'' is optional and should '...
               'be followed by its value if specified.']);
end

%% Find resonance peak of marker
B_target = 10E3 * g2b(marker_g,Pars.MWFQ);

% find 2 Gauss interval around suspected marker location
Interval = [B_target-5, B_target+5];
Slice = logical((Interval(1)<x) .* (x<Interval(2)));
x_slice = x(Slice);
y_slice = y(Slice);

if isempty(x_slice)
    disp('Marker location outside of spectrum. Proceeding without marker calibration.');
    xnew=x;
    return;
end

[~,~,~,p] = findpeaks(y_slice);
if p <0.5
    disp('No ESR signal found at marker location.');
    plot(x,y,x_slice,y_slice);
    xnew=x;
    return;
end

% determine apparent g-factor of marker
[~,B_marker] = gfactor_determination(x_slice,y_slice,Pars);

%% perform shift of x-axis
% calculate x-axis offset
B_offset = B_target - B_marker;
disp(['B_offset = ',num2str(B_offset),'G']);

% save new x-axis, confirm
xnew = x + B_offset;
disp('Offset corrected.');

% plot spectrum with marker location for visual confirmation
hold off; plot(xnew,y,'b');
hold on; plot(x,y,'b--');
plot(B_target,0,'o');plot(B_marker,0,'o');
plot(xnew,zeros(size(xnew)),'k','LineWidth',1);
axis tight;
legend('Corrected spectrum', 'Original spectrum','Corrected marker resonance','Original marker resonance');
hold off;
end