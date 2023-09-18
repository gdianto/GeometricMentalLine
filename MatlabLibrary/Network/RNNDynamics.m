function [t, Nu, Net] = RNNDynamics(Net, TimeSpan, NuExt, ODESolverOpts)

% Parameter processing.
if exist('NuExt','var') % External input.
   if ~isempty(NuExt)
      if size(Net.CParam.JExt,2) ~= size(NuExt,1) %JExt viene definito direttamente nel test
         error('Num. of rows in NuExt does not match the columns of Net.CParam.JExt')
      end
   end
else
   NuExt = [];
end

if numel(TimeSpan) == 1
   TimeSpan = [0 TimeSpan]+Net.t;
end

% Endogenous noise for each unit of the network.
if ~isfield(Net,'Noise')
   Net.Noise = [];
else
   if ~isempty(Net.Noise)
      if Net.Noise.SD > 0 % This condition is necessary to instantiate the endogenous noise.
         if ~isfield(Net.Noise,'dt') % Set default sampling step of noise.
            Net.Noise.dt = min(Net.MParam.Tau);
         end
         Net.Noise.t = TimeSpan(1):Net.Noise.dt:(TimeSpan(end)+Net.Noise.dt);
         % WARNING: This accelerates numerical integration but can overload the memory.
         if ~isfield(Net.Noise,'External') % Set endogenous or external noise.
            Net.Noise.External = true;
         end
         if Net.Noise.External
            Net.Noise.GWN = Net.Noise.SD * randn(numel(NuExt),numel(Net.Noise.t));
         else
            Net.Noise.GWN = Net.Noise.SD * randn(Net.P,numel(Net.Noise.t));
         end
      end
   end
end

OriginalIExt = Net.MParam.IExt;

if ~exist('ODESolverOpts','var')
   ODESolverOpts = [];
end

% Integration of the RNN dynamics.
[t,Nu] = ode45(@(t,Nu) odeRNN(t,Nu,Net,NuExt),TimeSpan,Net.Nu,ODESolverOpts);

% Reshape integration results.
t = t';
Nu = Nu';

% Restore original values and final state of the system...
Net.MParam.IExt = OriginalIExt;
Net.Nu = Nu(:,end);
Net.t = t(end);

end % [t, Nu, Net] = RNNDynamics(...)


%% Time gradient of the ODE to integrate the RNN.
%
function dNudt = odeRNN(t,Nu,Net,NuExt)

% Compute the contribution of the external input.
if ~isempty(NuExt)
   lNuExt = zeros(size(Net.CParam.JExt,2),1);
   if isstruct(NuExt)
      for k = 1:numel(lNuExt)
         lNuExt(k) = fnval(NuExt(k),t);
      end
   else
      for k = 1:numel(lNuExt)
         lNuExt(k) = NuExt{k}(t);
      end
   end
   Net.MParam.IExt = Net.CParam.JExt*lNuExt;
end

% Compute the additional endogenous and/or external noise.
if ~isempty(Net.Noise)
   ndx = ceil((t - Net.Noise.t(1))/Net.Noise.dt)+1;
   if Net.Noise.External
      Net.MParam.IExt = Net.MParam.IExt + Net.CParam.JExt*Net.Noise.GWN(:,ndx);
   else
      Net.MParam.IExt = Net.MParam.IExt + Net.Noise.GWN(:,ndx);
   end
end

% Compute the gradient.
dNudt = (Phi(Nu, Net) - Nu) ./ Net.MParam.Tau;

end
