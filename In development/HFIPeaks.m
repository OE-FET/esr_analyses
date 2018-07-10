function [locations, amplitudes]=HFIPeaks(I,A)
%% HFIPeaks
%
% Gives amplitudes and locations of peaks from isotropic HFI splitting in
% the limit of a strong magnetic field B >> A.
%
% INPUT:
% I - vector containing all nuclei with their spin
% A - corresponding HFI constants in magnetic field units
%
% OUTPUT:
% locations - peak locations around the resonance center (in units of A)
% amplitudes - peak amplitued in same order as locations
%
% EXAMPLE:
%
% HFIPeaks([1/2, 1],[3.6, 1.2])
% Computes resonance peaks for one nuclei with spin 1/2 and
% HFI constant of 3.6 and one nuclei with spin 1 and HFI constant 1.2 .
%
% WARNING: Computation time increases with 2^(sum(n)/20) .
%
%% 

nn=length(I); % total number of nuclei

for i=1:nn
    mN{i}=-I(i):I(i);
end

M=cell(nn,1);
[M{1:nn}]=ndgrid(mN{end:-1:1});
M=reshape(cat(nn+1,M{:}),[],nn);
M=M(:,end:-1:1); % all permutations of M_n values

Energy = M*A'; % Energies in units of HFI constant

% histogram of energies gives degeneracies and amplitudes
[amplitudes,locations]=hist(Energy,unique(Energy));

bar(locations,amplitudes,0.01);

end
