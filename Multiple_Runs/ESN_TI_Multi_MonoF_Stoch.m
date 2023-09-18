function [AnchorAcc, SymDist, MonkEffect, PCheck, Distr] = ESN_TI_Multi_MonoF_Stoch(Options, InputOptions)
% ESN in TI task

%input parameters
n_inputs = InputOptions.NInputs; % #images
waitT = InputOptions.Wait;
emissionT = InputOptions.Emission; % input signal width
emissionTR = InputOptions.EmissionR; % response signal width (+1=right,-1=left)
cycles = InputOptions.Cycles; % tMax ~ (wait+emission)*cycles

AT = 0;
AL = (waitT+emissionT)/Options.Tau * n_inputs*(n_inputs-1);
AP = (waitT+emissionT)/Options.Tau * n_inputs*(n_inputs-1); % (wait+emission) * (#images)*(#images-1)

AL = InputOptions.N_AL * AL;
SP = InputOptions.N_SP*AP;
AP = InputOptions.N_AP * AP;

AL = AL + 0.5*waitT/Options.Tau;

%-----discrete dynamics----
DISCRETE = Options.Discrete;
dt = Options.dt;
DELAY = InputOptions.Delay;
%-----
EFFECT_ONLY_OTH = true;

%-----continous dynamics----
INT_DENSITY = 10/Options.Tau; % n/Tau means integration_step=Tau/n
if(~DISCRETE)
    %DELAY = 100/INT_DENSITY % k/INT_DENS = k * integration_step
    DELAY = InputOptions.Delay;
end
%% Set the stimulation protocol.

OrthogonalStim = false;
%% Test protocol and other params.

TransientPeriod = AT*Options.Tau;
LearningPeriod = AL*Options.Tau;

TSpan = TransientPeriod + LearningPeriod;


NoiseForLearning = Options.NoiseForLearning; % size of the gaussian noise to make more 
                         % robust (attractive) the autonomous dynamics.
                         % It is the fraction of the std(Nu).
                         
                         
% Resampling allows to have better performances in computing the Wout.
% Indeed, adding a random reshuffling of the sampling time improves the
% quality of the pseudoinversion to obtain Wout.
ResampleLearnDyn = false;
%% Sets the network.

Net = setESN(Options);
%% Resampling

Net.t = 0;
Net.Nu = zeros(size(Net.Nu));

% Stimulation training + learning. Learning period is used to fit Wout.
% At the end, t = 0 is set to be the beginning of learning period.

% Settings to randomly resample the integrated dynamics: improves the fit of Wout.
if ResampleLearnDyn
    Dt = Options.Tau/5;
    TSpan = linspace(0,TSpan,ceil(TSpan/Dt)+1);
    TSpan(2:end-1) = TSpan(2:end-1) + (rand(1,numel(TSpan)-2)-0.5)*Dt/2; % (2,end-1) per mantanere gli estremi
end
%% Sets Inputs
epsT = 0.001 * Options.Tau;
interp_mode = 'linear';

% Only Adjacent couples of images
[InputAdj,ResponseAdj,~] = TI_Data_Multi_Mono(n_inputs,waitT,emissionT, emissionTR, epsT, cycles, interp_mode, 1);

% Free random couples
[Input,Response,tMax] = TI_Data_Multi_Mono_Regular(n_inputs,waitT,emissionT, emissionTR, epsT, cycles, interp_mode);

% singular image
[InputS,~] = TI_Data_Singular(n_inputs,waitT,emissionT, epsT, cycles, interp_mode);

%%
% smoothing data
smooth = true;
if(smooth)
    mn = 5 * (0.2 * emissionT/Options.Tau);
    x = linspace(0, tMax, 5*tMax/Options.Tau);
    
    mode = 'gaussian';
    for k=1:n_inputs
        
        InputAdj{k} = griddedInterpolant(x,...
            smoothdata(InputAdj{k} (x), mode, mn),interp_mode);
        
        
        Input{k} = griddedInterpolant(x,...
            smoothdata(Input{k} (x), mode, mn),interp_mode);
        
        
        InputS{k} = griddedInterpolant(x,...
            smoothdata(InputS{k} (x), mode, mn),interp_mode);
        

    end
    
    ResponseAdj{1} = [griddedInterpolant(x,...
    smoothdata(ResponseAdj{1} (x), mode, mn),interp_mode)];

    Response{1} = [griddedInterpolant(x,...
    smoothdata(Response{1} (x), mode, mn),interp_mode)];
end
%%
% Shift Response
ResponseAdj{1} = @(t) ResponseAdj{1} (t-DELAY);
Response{1} = @(t) Response{1} (t-DELAY);
%%
NuExtAdj = [InputAdj, ResponseAdj]';
NuExt = [Input, Response]';
NuExtS = [InputS, Response]';
%%
NuExtTrans = NuExtAdj;
NuExtLearn = NuExtAdj;
NuExtAuto = NuExt;
%% Sets the stimulation protocol.

if OrthogonalStim
    Net.CParam.JExt = 2*(rand(Net.P,size(NuExtLearn, 1)) > 0.5) - 1; % Orthogonal (discrete spin) stimulation. (+-1)
    Net.CParam.JExt(:,1:n_inputs) = Net.CParam.JExt(:,1:n_inputs)*Options.StdJinput;
    Net.CParam.JExt(:,n_inputs+1:end) = Net.CParam.JExt(:,n_inputs+1:end)*Options.StdJback;
else
    Net.CParam.JExt = randn(Net.P,size(NuExtLearn, 1));
    Net.CParam.JExt(:,1:n_inputs) = Net.CParam.JExt(:,1:n_inputs)*Options.StdJinput;
    Net.CParam.JExt(:,n_inputs+1:end) = Net.CParam.JExt(:,n_inputs+1:end)*Options.StdJback;
end
%% add Directional Stimulous
Net.DirNoise.Intensity = Options.DirNoise.Intensity;
Net.DirNoise.ENdx = n_inputs+1;
%% Numerically integrate the network dynamics
% Transiente + Apprendimento

opts = odeset('MaxStep',(0.5*Options.Tau));

if(ResampleLearnDyn)
    IP_indx = find(TSpan <= TransientPeriod);
    LP_indx = find(TSpan > TransientPeriod);
    if(DISCRETE)
        [~, ~, Net] = RNNDynamicsDisc(Net, TSpan(IP_indx), NuExtTrans, dt);
        Net.t=0;
        [tLP, NuLP, Net] = RNNDynamicsDisc(Net, TSpan(LP_indx), NuExtLearn, dt);
    else
        [~, ~, Net] = RNNDynamics(Net, TSpan(IP_indx), NuExtTrans, opts);
        Net.t=0;
        [tLP, NuLP, Net] = RNNDynamics(Net, TSpan(LP_indx), NuExtLearn, opts);
    end
    
else
    TTL = linspace(0, TransientPeriod, ceil(INT_DENSITY*TransientPeriod));
    TLL = linspace(0, LearningPeriod, ceil(INT_DENSITY*LearningPeriod));
    if(DISCRETE)
        [~, ~, Net] = RNNDynamicsDisc(Net, TransientPeriod, NuExtTrans, dt);
        Net.t=0;
        [tLP, NuLP, Net] = RNNDynamicsDisc(Net, LearningPeriod, NuExtLearn, dt);
    else
        [~, ~, Net] = RNNDynamics(Net, TTL, NuExtTrans, opts);
        Net.t=0;
        [tLP, NuLP, Net] = RNNDynamics(Net, TLL, NuExtLearn, opts);
    end
end
%% Learn the coefficients "Wout"

% This is without resampling the integrated dynamics.
OutOrig = NuExtLearn{1+n_inputs}(tLP);

if NoiseForLearning > 0
   AddedNoise = NoiseForLearning*std(NuLP(:))*randn(size(NuLP));
   %Wout = OutOrig*pinv(NuLP'+AddedNoise')';
   Wout = OutOrig/(NuLP+AddedNoise);
else
   %Wout = OutOrig*pinv(NuLP')';
   Wout = OutOrig/NuLP;
end
%% Incorporates the learnt synaptic matrix and test the autonomous dynamics.

NetSelf = Net;
NetSelf.CParam.JOpen = NetSelf.CParam.J;
NetSelf.CParam.JToClose = Net.CParam.JExt(:,n_inputs+1:end)*Wout;

if(Options.StdJback > 0)
    NetSelf.CParam.J = NetSelf.CParam.JOpen + NetSelf.CParam.JToClose;
    NetSelf.CParam.JExt(:,n_inputs+1:end) = 0 * NetSelf.CParam.JExt(:,n_inputs+1:end);
else
    NetSelf.CParam.J = NetSelf.CParam.JOpen;
end

AutonomousPeriod = AP*Options.Tau;
TAA = linspace(LearningPeriod, LearningPeriod+AutonomousPeriod, ceil(INT_DENSITY*AutonomousPeriod));
if(DISCRETE)
    [tA, NuA, NetSelf] = RNNDynamicsDisc(NetSelf, AutonomousPeriod, NuExtAuto, dt);  
else
    [tA, NuA, NetSelf] = RNNDynamics(NetSelf, TAA, NuExtAuto, opts); 
end
%%
% Read out.
OutA = Wout*NuA;
OutAOrig(1,:) = NuExtAuto{1+n_inputs}(tA);

%% Event Definition
event = [];
counter = 0;
Sdist.counter = zeros(1,n_inputs-1); % num. events for each symb. distance
lastResp = 0;
eventProp.sdist = 0;
tStart = tA(1);
minOutMax = inf;

relaxingCond = zeros(size(tA));
for t=tA
    if((NuExtAuto{end} (t) ~= 0) && (lastResp == 0)) % event starts
        tStart = t - DELAY;
        counter = counter + 1;
        event = [event {counter}];
        
        sd = 0;
        for k=1:n_inputs
            boolDx = ceil(0.9*NuExtAuto{k} (t));
            boolSx = ceil(-0.9*NuExtAuto{k} (t));
            sd = sd + k*boolDx - k*boolSx;
            
            if(boolDx)
                eventProp.image.dx = k;
            elseif(boolSx)
                eventProp.image.sx = k;
            end
            
        end
        
        eventProp.sdist = abs(sd);
        Sdist.counter(eventProp.sdist) = Sdist.counter(eventProp.sdist) +1;
        
        if(t > waitT)
            relaxingCond = relaxingCond + ((tA < tStart) & (tA > tStart - 0.5*waitT));
        end
    elseif( ((NuExtAuto{end} (t) == 0) && (lastResp ~= 0)) || ... % event ends
            ((t == tA(end)) && (lastResp ~= 0)) )
        interval = find((tA >= tStart) & (tA <= t));
        eventProp.t = [tStart t];
        eventProp.Out = OutA(interval);
        [~,maxNdx] = max(abs(eventProp.Out));
        eventProp.maxOut = eventProp.Out(maxNdx);
        
        if(abs(eventProp.maxOut) < minOutMax)
            minOutMax = abs(eventProp.maxOut);
        end
        
        event{counter} = eventProp;
    end
    lastResp = NuExtAuto{end} (t);
    
end
%% Accuracy - MaxOutput - Reaction Times - Anchor Effect
th = 0.;
if(Options.th == 999)
    th = 0.95*minOutMax;
else
    th = Options.th;
end

%fprintf("TH="+num2str(th) +"\n")

Sdist.acc = zeros(1, n_inputs-1);
Sdist.maxOut = zeros(1, n_inputs-1);
Sdist.reactTimes = zeros(1, n_inputs-1);

AccHeat = zeros(n_inputs,n_inputs);
counterHeat = zeros(n_inputs,n_inputs);
counter_oTH_Heat = zeros(n_inputs,n_inputs);
SD_oTH_counter = zeros(1,n_inputs-1);
for k=1:numel(event)
    % maximum output
    Sdist.maxOut(event{k}.sdist) = Sdist.maxOut(event{k}.sdist) + abs(event{k}.maxOut);

    % reaction times
    interval = find((tA >= event{k}.t(1)) & (tA <= event{k}.t(end)));
    tInt = tA(interval);
    OutInt = OutA(interval);
    OutOrigInt = OutAOrig(interval);
    
    guess = 0;
    if(~isempty(find(abs(OutInt) > th)))
        ndxOTH = find(abs(OutInt) > th);
        tTH = tInt(ndxOTH);
        Sdist.reactTimes(event{k}.sdist) = Sdist.reactTimes(event{k}.sdist) + (tTH(1)-event{k}.t(1)); 
    
        %accuracy
        guess = (sign(OutInt(ndxOTH(1))*OutOrigInt(ndxOTH(1))) > 0);

        counter_oTH_Heat(event{k}.image.dx,event{k}.image.sx) = counter_oTH_Heat(event{k}.image.dx,event{k}.image.sx) + 1;
        SD_oTH_counter(event{k}.sdist) = SD_oTH_counter(event{k}.sdist) + 1; %just for reaction times
    end
    Sdist.acc(event{k}.sdist) = Sdist.acc(event{k}.sdist) + guess;

    % accuracy heatmap
    AccHeat(event{k}.image.dx,event{k}.image.sx) = AccHeat(event{k}.image.dx,event{k}.image.sx) + guess;
    counterHeat(event{k}.image.dx,event{k}.image.sx) = counterHeat(event{k}.image.dx,event{k}.image.sx) + 1;

end


SymDist.Acc = zeros(1, n_inputs-1);
%SymDist.Acc = Sdist.acc ./ Sdist.counter;
SymDist.MaxOut = Sdist.maxOut ./ Sdist.counter;
SymDist.ReactTimes = Sdist.reactTimes ./ SD_oTH_counter;
SymDist.Counter = Sdist.counter;
SymDist.OTHCounter = SD_oTH_counter;

if(EFFECT_ONLY_OTH)
    SymDist.Acc = Sdist.acc ./ SD_oTH_counter;
    AccHeat = AccHeat ./ counter_oTH_Heat;
else
    SymDist.Acc = SymDist.Acc ./ Sdist.counter;
    AccHeat = AccHeat ./ counterHeat;
end


AnchorAcc = zeros(1,n_inputs);
AccHeatRestricted = AccHeat(2:end-1,2:end-1);
for k=1:n_inputs
    if(k==1 || k==n_inputs)
        nC = 2*numel(AccHeat(k,~isnan(AccHeat(k,:))));
        AnchorAcc(k) = (sum(AccHeat(k,~isnan(AccHeat(k,:)))) + sum(AccHeat(~isnan(AccHeat(:,k)), k)))/nC;
    else
        nC = 2*(numel(AccHeatRestricted(k-1,~isnan(AccHeatRestricted(k-1,:)))));
        AnchorAcc(k) = (sum(AccHeatRestricted(k-1,~isnan(AccHeatRestricted(k-1,:)))) + sum(AccHeatRestricted(~isnan(AccHeatRestricted(:,k-1)),k-1)))/nC;
    end
end
%% Correlation
sdist = zeros(size(OutAOrig));
for k=1:n_inputs
    sdist = sdist + k*ceil(0.9*NuExtAuto{k} (tA)) - k*ceil(-0.9*NuExtAuto{k} (tA));
end
sdist = abs(sdist);

%corrTh = 0.4;
neuroPercs = zeros(1,n_inputs-1);
for k=1:n_inputs-1
    M = [OutA(sdist==k);NuA(:,sdist==k)]';
    Rm = corrcoef(M);
    Rc = Rm(1,2:end);
    neuroPercs(k) = mean(abs(Rc));
end
SymDist.NeuroPercs = neuroPercs;
%% Controllo Percettrone
% evar pcs
[~,~,EVs_A] = pca(NuA');
eVar = cumsum(EVs_A)/sum(EVs_A)*100;
PCheck.EV_PC1 = eVar(1);
%%
% Linearità Wout su Jext
J = NetSelf.CParam.JExt(:,1:end-1);
pj = Wout * J;
pj = (pj - pj(1));
pj = pj / max(abs(pj)) * (n_inputs-1);

x = 1:n_inputs;
c = polyfit(x,pj,1);
y_est = polyval(c,x);

err_pj_temp = (y_est - pj).^2;
PCheck.Err_Pj = sqrt(sum(err_pj_temp))/n_inputs;
%% Effects Controll

monkAccBool = true;
for sd = 1:(n_inputs-2)
    monkAccBool = monkAccBool && (SymDist.Acc(sd) < SymDist.Acc(sd+1));
end
MonkEffect.AccBool = monkAccBool;

monkRtBool = true;
for sd = 1:(n_inputs-2)
    monkRtBool = monkRtBool && (SymDist.ReactTimes(sd) > SymDist.ReactTimes(sd+1));
end
MonkEffect.RtBool = monkRtBool;

monkAncBool = true;
for sd = 1:n_inputs-2
    monkAncBool = monkAncBool && (AnchorAcc(1) > AnchorAcc(sd+1));
    monkAncBool = monkAncBool && (AnchorAcc(end) > AnchorAcc(end-sd));
end
MonkEffect.AncBool = monkAncBool;

monkCorrBool = true;
for sd = 1:(n_inputs-2)
    monkCorrBool = monkCorrBool && (SymDist.NeuroPercs(sd) < SymDist.NeuroPercs(sd+1));
end
MonkEffect.CorrBool = monkCorrBool;

MonkEffect.All = (monkAccBool & monkRtBool & monkAncBool & monkCorrBool);

%% Singular Stimulation

SingPeriod = SP*Options.Tau;

TPP = linspace(AutonomousPeriod, AutonomousPeriod+SingPeriod, ceil(INT_DENSITY*SingPeriod));
if(DISCRETE)
    [tS, NuS, ~] = RNNDynamicsDisc(NetSelf, SingPeriod, NuExtS, dt);  
else
    [tS, NuS, ~] = RNNDynamics(NetSelf, TPP, NuExtS, opts);
end

OutRank = Wout*NuS;
% Event Singular Definition

event = [];
counter = 0;
lastResp = 0;
tStart = tS(1);
minOutMax = inf;

relaxingCond = zeros(size(tS));
for t=tS
    if((NuExtS{end} (t) ~= 0) && (lastResp == 0)) % event starts
        tStart = t - DELAY;
        counter = counter + 1;
        event = [event {counter}];
        
        for k=1:n_inputs
            boolDx = ceil(0.9*NuExtS{k} (t));
            
            if(boolDx)
                eventProp.image.dx = k;
            end
        end
        
        
        if(t > waitT)
            relaxingCond = relaxingCond + ((tS < tStart) & (tS > tStart - 0.5*waitT));
        end
    elseif( ((NuExtS{end} (t) == 0) && (lastResp ~= 0)) || ... % event ends
            ((t == tS(end)) && (lastResp ~= 0)) )
        eventProp.t = [tStart t];
        eventProp.Out = OutRank((tS >= tStart) & (tS <= t));
        [~,maxNdx] = max(abs(eventProp.Out));
        eventProp.maxOut = eventProp.Out(maxNdx);
        
        if(abs(eventProp.maxOut) < minOutMax)
            minOutMax = abs(eventProp.maxOut);
        end
        
        event{counter} = eventProp;
    end
    lastResp = NuExtS{end} (t);
end

%%
distrImg = nan(n_inputs, numel(event));
for k=1:numel(event)
    distrImg(event{k}.image.dx,k) = event{k}.maxOut;
end
%%
RM = round(numel(event{1}.Out)/2.3);
distrImgContinous = nan(n_inputs, numel(event)*numel(event{1}.Out));
for k=1:numel(event)
    centralSignal = event{k}.Out(RM:end-RM);
    distrImgContinous(event{k}.image.dx,((k-1)*numel(centralSignal)+1):k*numel(centralSignal)) = centralSignal;
end
%% Wilcoxon Test
wrs_pv = nan(n_inputs);
wrsC_pv = nan(n_inputs);
for k=1:n_inputs
    for v=1:n_inputs
        if(k~=v)
            A = distrImg(k,~isnan(distrImg(k,:)));
            B = distrImg(v,~isnan(distrImg(v,:)));
            wrs_pv(k,v) = ranksum(A,B);
            A = distrImgContinous(k,~isnan(distrImgContinous(k,:)));
            B = distrImgContinous(v,~isnan(distrImgContinous(v,:)));
            wrsC_pv(k,v) = ranksum(A,B);
        end
    end
end
%%
wrs_im = zeros(1,n_inputs);
for k=1:n_inputs
    nC = 2*numel(wrs_pv(k,~isnan(wrs_pv(k,:))));
    wrs_im(k) = (sum(wrs_pv(k,~isnan(wrs_pv(k,:)))) + sum(wrs_pv(~isnan(wrs_pv(:,k)), k)))/nC;
end

wrsC_im = zeros(1,n_inputs);
for k=1:n_inputs
    nC = 2*numel(wrsC_pv(k,~isnan(wrsC_pv(k,:))));
    wrsC_im(k) = (sum(wrsC_pv(k,~isnan(wrsC_pv(k,:)))) + sum(wrsC_pv(~isnan(wrsC_pv(:,k)), k)))/nC;
end
%%
Distr.wrs_im = wrs_im;
Distr.wrsC_im = wrsC_im;