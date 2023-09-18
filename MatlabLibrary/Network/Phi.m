function NuOut = Phi(NuIn, Net)

% sigmoid = @(x)1./(1+exp(-4*x)); % Logistic function as a shifted hyperbolic tangent.
sigmoid = @(x)tanh(x);            % As in (Sompolisky et al., PRL 1988).

if size(NuIn,1) == 1
   NuIn = NuIn';
end

NuOut = sigmoid(Net.CParam.J*NuIn + Net.MParam.IExt)./Net.MParam.Tarp;
