function CM = gradedColormap(StartingColor,EndingColor,CMHalfSize)
%
% CM = gradedColormap(StartingColor,EndingColor[,CMHalfSize])
%
%   Copyright Maurizio Mattia - Mar. 14, 2021
%

if ~exist('CMHalfSize','var')
   CMHalfSize = 32;
end
FirstCM = linspace(1,0,CMHalfSize)'*StartingColor;
FirstCM = FirstCM + linspace(0,1,CMHalfSize)'*[1 1 1];
LastCM = linspace(0,1,CMHalfSize)'*EndingColor;
LastCM = LastCM + linspace(1,0,CMHalfSize)'*[1 1 1];
CM = [FirstCM; LastCM(2:end,:)];
