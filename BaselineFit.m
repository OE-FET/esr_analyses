function [ycorr, yfit] = BaselineFit(y)
%BASELINEFIT Performs a baseline fit on the input data
%   [ycorr, yfit] = BASELINEFIT(y) performs a baseline fit on the input data.
%   The baseline regions can be be selected through a GUI, and the baseline is
%   fitted as a spline through the smoothed data. Adjust the number points 
%   to smooth over according to the noise in the data.
%
%   OUTPUT:
%   ycorr - baseline corrected spectrum
%   yfit - fitted baseline
%
%   DEPENDENCIES:
%   StackPlot.m
%   PointInput.m
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2018/07/05 12:58 $    $Revision: 1.1 $

%% process data
% get spectrum size
dim = size(y);
lst = dim(1);
newdimy = [dim(1), prod(dim(2:end))];
y = reshape(y, newdimy);
x = 1:lst;

%% select methods for interpolation
% use spline interpolation for polynomial baseline fit
method  = 'spline';
% average over 10 points for smoothing before fit 
avgpts  = round(lst/100); % 1/100 of length of data 
step = 8*avgpts;

%% select area for baseline fit

ok = false;
% re-promt for user input until user confirms good baseline correction
while ~ok
    % plot spectrum
    bffig = StackPlot(x, real(y));

    % set title of plot
    title('Baseline Fit - Select area of spectrum');
    % promt user for input of baseline areas
    fprintf(['\n Now select the area of the spectrum,',...
        '\n by indicating points with the courser.',...
        '\n Press Enter key when done.\n'])
    % open GUI for input, accept only 2 points
    [a, b] = PointInput(2);                                 %#ok
    % round to integer values
    bounds = round(a.');
    bounds = sort(bounds);
    
    % if points outside x-axis are selected, replace them with axis limits
    if bounds(1)<1,       bounds(1)=1;            end     
    if bounds(end)>lst,   bounds(end)=lst;        end     
    
    % do nothing if whole x-axis range is selected
    if (bounds(1)<1 &&bounds(end)>lst)
        ycorr = y; yfit = zeros(size(y));
        return;
    end 
    
    % within baseline areas, use points in intervals of step for polynomial
    % fitting
    pts = [(1:step:bounds(1)), (bounds(2):step:lst)];
    pts = round(pts);
    npts = numel(pts);
    
    % smooth curve by averaging over avgpts neighboring points
    pss = zeros(npts, 2);
    pss(:,1) = pts - floor(avgpts/2);
    pss(:,2) = pss(:,1) + avgpts;
    pss(pss < 1) = 1;
    pss(pss > lst) = lst;
    yavg = zeros([npts, newdimy(2)]);
    for n = 1:npts
        yavg(n,:) = mean(y(pss(n,1):pss(n,2),:),1);
    end
    %% perform baseline fit
    yfit = interp1(pts, yavg, x, method);
    if size(yfit, 1)==1
        yfit = shiftdim(yfit, 1);    % make yfit a column if it is a row vector
    end
%     yfit(x<bounds(1), :) = y(x<bounds(1), :);
%     yfit(x>bounds(4), :) = y(x>bounds(4), :);
    
    % plot for visual confirmation
    hold on;
    yL = get(gca,'YLim');
    [~, phandle] = StackPlot(x, real(yfit), 'yoffset', 0.5*max(max(y)));    
    set(phandle, 'Color', 'blue');
    ylim(gca, yL);
    for i = 1:2
        line([bounds(i) bounds(i)], yL, 'Color', 'c');
    end
    % plot(x,real(y),'b',x,real(yfit),'r');
    hold off;

    %% promt user for confirmation of fit
    set(bffig,'Name','Baseline Fit - Verify baseline')
    answer = input('  Do you to redo fit and reselect baseline points?[N] ','s');
    if isempty(answer),     answer = 'n';   end
    if strcmpi(answer,'y')
        ok = false;
    else
        ok = true;
    end
end
% close figure if it exists
if any(findobj('Type','figure')==bffig)
    close(bffig),                   
end
% subtact baseline 
ycorr = y - yfit;
% reshape spectrum to original form
ycorr = reshape(ycorr,dim);
yfit = reshape(yfit,dim);