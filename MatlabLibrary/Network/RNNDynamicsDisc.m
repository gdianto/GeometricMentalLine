function [TimeSpan, NuArr, Net] = RNNDynamicsDisc(Net, TimeSpan, NuExt, dt)
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
   TimeSpan = Net.t:dt:Net.t+TimeSpan;
end

% Endogenous noise for each unit of the network.
if ~isfield(Net,'Noise')
   Net.Noise = [];
else
   if ~isempty(Net.Noise)
      if Net.Noise.SD > 0 % This condition is necessary to instantiate the endogenous noise.
         if ~isfield(Net.Noise,'dt') % Set default sampling step of noise.
            Net.Noise.dt = dt;
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

%---
if ~isfield(Net, 'DirNoise')
    Net.DirNoise.Intensity = [];
    Net.DirNoise.ENdx = [];
else
    if ~isfield(Net.DirNoise, 'Intensity')
        Net.DirNoise.Intensity = [];
    end
    
    if ~isfield(Net.DirNoise, 'ENdx')
        Net.DirNoise.ENdx = [];
    end
end
%---
        

OriginalIExt = Net.MParam.IExt;

% Integration of the RNN dynamics.
NuArr = DiscDyn(TimeSpan, Net, NuExt);

% Reshape integration results.
%t = t';
%NuArr = NuArr';

% Restore original values and final state of the system...
Net.MParam.IExt = OriginalIExt;
Net.Nu = NuArr(:,end);
Net.t = TimeSpan(end);

end % [t, Nu, Net] = RNNDynamics(...)


function NuArr = DiscDyn(allTimes,Net,NuExt)

NuArr = nan(Net.P, numel(allTimes));
Nu = Net.Nu;
NuArr(:,1) = Nu;

%-----
if ~isempty(Net.DirNoise.Intensity)
    dJExt = randn(size(Net.CParam.JExt)) * Net.DirNoise.Intensity;
    if ~isempty(Net.DirNoise.ENdx)
        dJExt(:,Net.DirNoise.ENdx) = 0 * dJExt(:,Net.DirNoise.ENdx);
    end
    JExt_N = Net.CParam.JExt + dJExt;
    normJE = (std(Net.CParam.JExt) + 0.001) ./ (std(JExt_N) + 0.001);
    JExt_N = JExt_N .* normJE;
    changed = false;
end
%------

    for v= 2 : numel(allTimes)
        t = allTimes(1,v);
        dt = t - allTimes(1,v-1);
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
           
           
           % directional noise
            if ~isempty(Net.DirNoise.Intensity)
                if(Net.DirNoise.Intensity > 0)
                    if(max(abs(Net.MParam.IExt)))
                        Net.MParam.IExt = JExt_N * lNuExt;
                        changed = false;
                    else
                        if(~changed)
                            
                            if isfield(Net,'Patterns')
                                dDirCos = randn(size(Net.Patterns,2),size(NuExt, 1));
                                dDirCos = dDirCos/norm(dDirCos);
                                dJExt = Net.Patterns * dDirCos * Net.DirNoise.Intensity;
                            else
                                dJExt = randn(size(Net.CParam.JExt)) * Net.DirNoise.Intensity;
                            end
                            
                            
                            if ~isempty(Net.DirNoise.ENdx)
                                dJExt(:,Net.DirNoise.ENdx) = 0 * dJExt(:,Net.DirNoise.ENdx);
                            end
                            JExt_N = Net.CParam.JExt + dJExt;
                            normJE = (std(Net.CParam.JExt) + 0.001) ./ (std(JExt_N) + 0.001);
                            JExt_N = JExt_N .* normJE;
                            changed = true;
                        end
                    end
                end
            end
            %--
            
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
        Nu = dt/Net.MParam.Tau(1) * (Phi(Nu, Net)-NuArr(:,v-1)) + NuArr(:,v-1);
        NuArr(:,v) = Nu;
    end

end
