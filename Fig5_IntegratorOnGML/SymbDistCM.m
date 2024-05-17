function CM = SymbDistCM(NoSD)
%
%   Colormap used to represent differen Symbolic Distance in  
%   (Di Antonio et al., bioRxiv 2023, DOI: 10.1101/2023.08.03.551859).
%
%   Copyright Maurizio Mattia - Mar. 3, 2023
%

PosSD = [255 80 155]/255;
NegSD = [70 180 255]/255;
% NoSD = 4;

CM = flipud((PosSD' * linspace(1,0,2*NoSD))' + (NegSD' * linspace(0,1,2*NoSD))');