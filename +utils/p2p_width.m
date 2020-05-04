function width = p2p_width(x, y)
% P2P_WIDTH returns the peak-to-peak width of a derivative resonance line
%   in the same units as x.

[~, i0] = max(y);
[~, i1] = min(y);

width = x(i1) - x(i0);

end