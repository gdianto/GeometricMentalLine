function [InputF,ResponseF,tMax] = TI_Data_Multi_Mono_Regular(n_inputs,waitT,emissionT, emissionTR, epsT, cycles, interpMode)
%   TI data Generator:
%   Input is a matrix whose raws represent the input channels
%


% time sequence over which interpolate data
deltaT = [epsT, emissionT, epsT, waitT];
deltaT = repmat(deltaT, 1, cycles);
deltaT = [0, waitT, deltaT];
t = cumsum(deltaT);
tMax = t(end);

% time sequence over which interpolate data
waitTR = waitT - (emissionTR - emissionT);
deltaTR = [epsT, emissionTR, epsT, waitTR];
deltaTR = repmat(deltaTR, 1, cycles);
deltaTR = [0, waitTR, deltaTR];
tR = cumsum(deltaTR);


% input matrix generator
Input = zeros(n_inputs, numel(t));
Response = zeros(1, numel(t));

counter = 0;
counterLim = n_inputs;
imgList = 1:n_inputs;
r = 1;
for k=4:4:numel(t) % every 4 index we find an emission (between k-1 and k)
    
    counter = counter +1;
    if(counter == r)
        counter = r+1;
    end
    
    if(counter == n_inputs+1)
        counter = 1;
        r = r +1;
        
        if(r == n_inputs+1) % restart the cycle
            counter = 2;
            r = 1;
        end
        
    end
    r2 = imgList(counter);
    
    Input(r,k) = 1;
    Input(r2,k) = -Input(r,k);
    Input(r,k-1) = Input(r,k);
    Input(r2,k-1) = Input(r2,k);
    
    if(r>r2)
        rH = r;
    else
        rH = r2;
    end
    
    Response(1,k) = Input(rH,k);
    Response(1,k-1) = Response(1,k);
end

% interpolate
InputF = [];
ResponseF = [];
for k=1:n_inputs
    InputF = [InputF {k}];
    InputF{k} = griddedInterpolant(t,Input(k,:),interpMode);
end
ResponseF = [ResponseF {1}];
ResponseF{1} = griddedInterpolant(tR,Response(1,:), interpMode);

end

