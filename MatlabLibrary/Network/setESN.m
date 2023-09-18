function Net = setESN(Options)
% The recurrent neural network (RNN) we are setting is r-based rather than v-base 
% referring to the notation introduce in (Miller & Fumarola, NC 2012). 
% Although these two dynamics have been proven to be equivalent (Miller & 
% Fumarola, NC 2012), the rate (r-)based is the one introduced in (Wilson & 
% Cowan, BJ 1972) and is the closest to the low-D mean-field dynamics of 
% spiking neuron networks (Mattia, Arxiv 2016).
% As single-unit gain function we use tanh(), meaning that firing rates
% range between -1.0 and 1.0.


%% Network parameters.
%
NetSize = Options.NetSize;       % Number of units/modules (N).

if ~isfield(Options,'Tau') % Single-unit decay constant (inertia) in seconds.
   Options.Tau = 1.0; 
end
Tau = Options.Tau;

if ~isfield(Options,'MaxNu') % This is the maximum output rate and multiplies Phi().
   Options.MaxNu = 1.0; 
end
MaxNu = Options.MaxNu;

if ~isfield(Options,'UnitSize') % Neurons per module (magnitude of finite-size noise).
   Options.UnitSize = Inf; 
end
UnitSize = Options.UnitSize;

if ~isfield(Options,'Iext') % Offset of the input field, a kind of external input current summed to the recurrent field.
   Options.Iext = 0.0; 
end
Iext = Options.Iext;

if ~isfield(Options,'MeanJself') % 's' param of (Stern et al, 2014), s > 1 -> bistable module.
   Options.MeanJself = 0.0; 
end
MeanJself = Options.MeanJself;

if ~isfield(Options,'StdJinter') % 'g' param of (Stern et al, 2014) and (Sompolinsky et al., PRL 1988).
   Options.StdJinter = 1.0;      % The paramagnetic state is marginally stable.
end
StdJinter = Options.StdJinter;

if ~isfield(Options,'EqualizeJinter') % if 1, sets the rows of the synaptic matrix to have 0 mean and g as standard deviation.
   Options.EqualizeJinter = 0;      
end
EqualizeJinter = Options.EqualizeJinter;

if ~isfield(Options,'NoiseSD') % Are endogenous noise params provided?
   Options.NoiseSD = 0;
end

%% Network setting with random definition ...
%
Net.P = NetSize; % Number of units/modules.
Net.MParam.N    = repmat(UnitSize,Net.P,1);
Net.MParam.Tau  = repmat(Tau,Net.P,1); % In ms.
Net.MParam.Tarp = repmat(1/MaxNu,Net.P,1);   % 1/Tarp set the maximum firing rate.
Net.MParam.IExt = repmat(Iext,Net.P,1);          % External current.


Net.CParam.J = randn(Net.P,Net.P);
if EqualizeJinter
   for m = 1:Net.P
      ndx = [1:m-1 m+1:Net.P]; %ovvero prendo tutti tranne m (evito la diagonale)
      Net.CParam.J(m,ndx) = Net.CParam.J(m,ndx) - mean(Net.CParam.J(m,ndx));
      Net.CParam.J(m,ndx) = Net.CParam.J(m,ndx)/std(Net.CParam.J(m,ndx));
   end
end

Net.CParam.J = Net.CParam.J .* sqrt(repmat(Net.MParam.Tarp/Net.P,1,Net.P)); % normalization
Net.CParam.J = Net.CParam.J*StdJinter/max(abs(eig(Net.CParam.J)));


ndxSelf = 1:Net.P+1:Net.P^2; %1 1+(1+P) 1+(1+P)+(1+P) ... P^2 , sarebbero gli indici della diagonale
Net.CParam.J(ndxSelf) = ones(size(ndxSelf))*MeanJself;
Net.CParam.J(ndxSelf) = Net.CParam.J(ndxSelf).*Net.MParam.Tarp';

Net.StdJinter = StdJinter;
Net.MeanJself = MeanJself;

%% Set the parameters of the endogenous noise, if any.
%
if Options.NoiseSD > 0
   Net.Noise.SD = Options.NoiseSD;
   if isfield(Options,'NoiseDt')
      Net.Noise.dt = Options.NoiseDt;
   end
else
   Net.Noise = []; % No noise.
end

%% Set the initial condition of the network dynamics (t = 0 and Nu = 0).
%
Net.t = 0;
Net.Nu = zeros(Net.P,1);
