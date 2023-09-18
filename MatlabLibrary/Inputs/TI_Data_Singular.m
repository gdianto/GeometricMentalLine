function [InputF,tMax] = TI_Data_Singular(n_inputs,waitT,emissionT, epsT, cycles, interpMode)
%   TI data Generator:
%   Input is a matrix whose raws represent the input channels
%


% time sequence over which interpolate data
deltaT = [epsT, emissionT, epsT, waitT];
deltaT = repmat(deltaT, 1, cycles);
deltaT = [0, waitT, deltaT];
t = cumsum(deltaT);
tMax = t(end);




% input matrix generator
Input = zeros(n_inputs, numel(t));

counter = 0;
imgList = 1:n_inputs;
r = 1;
for k=4:4:numel(t) % every 4 index we find an emission (between k-1 and k)
    Input(r,k) = 1;
    Input(r,k-1) = Input(r,k);
    r = r+1;
    if(r == n_inputs+1)
        r = 1;
    end
end

% interpolate
InputF = [];
for k=1:n_inputs
    InputF = [InputF {k}];
    InputF{k} = griddedInterpolant(t,Input(k,:),interpMode);
end

end

