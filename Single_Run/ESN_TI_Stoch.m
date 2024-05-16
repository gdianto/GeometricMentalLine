%% TI task
clear all
close all
%%
tic
seed = randi(9999);
rng(seed);

PlotFigures = true;

% Set the network.
Options.NetSize = 100;   % Number of units.
Options.Tau = 0.1;     % Single-unit decay constant in seconds.
Options.StdJinter = 0.05; % 'g' param of (Stern et al, 2014) and (Sompolinsky et al., PRL 1988).
OrthogonalStim = false;
Options.StdJinput = Options.StdJinter*0.2;
Options.StdJback = Options.StdJinput*2; %1.8

Options.NoiseSD = 0.0;% Standard deviation of the Gaussian white noise added to Iext.
% Options.NoiseDt = Options.Tau/2;   % Sampling time of endogenous noise (default Tau).
NOISE_INT = 0.0; % Directional Noise Intensity
Options.th = 0.5; % response threshold

EFFECT_ONLY_OTH = true;

%input parameters
N_INPUTS = 7; % #images
WAIT = 30*Options.Tau; % waiting time between trials
EMISSION = 15*Options.Tau; % input signal width
EMISSION_R = 15*Options.Tau; % response signal width (+1=right,-1=left)
CYCLES = 10000; % tMax ~ (wait+emission)*cycles

AT = 0;
AL = (WAIT+EMISSION)/Options.Tau * N_INPUTS*(N_INPUTS-1);
AP = (WAIT+EMISSION)/Options.Tau * N_INPUTS*(N_INPUTS-1); % (wait+emission) * (#images)*(#images-1)

N_AL = 20 / N_INPUTS;
AL = N_AL * AL + 0.5*WAIT/Options.Tau;

N_SP = 30 / N_INPUTS; %30
SP = N_SP * AP;

N_AP = 20;
AP = N_AP * AP;

%-----discrete dynamics options----
dt = 0.1*Options.Tau;
DELAY = 1*dt;
%-----

NoiseForLearning = 0.05;
%% Test protocol and other params.

TransientPeriod = AT*Options.Tau;
LearningPeriod = AL*Options.Tau;

TSpan = TransientPeriod + LearningPeriod;
%% Sets the network.

Net = setESN(Options);
%% Sets Inputs

n_inputs = N_INPUTS;
waitT = WAIT;
emissionT = EMISSION;
emissionTR = EMISSION_R;
epsT = 0.001 * Options.Tau;
cycles = CYCLES; % tMax ~ (wait+emission)*cycles
interp_mode = 'linear';

% Only adjacent couples of images
[InputAdj,ResponseAdj,tMaxAdj] = TI_Data_Multi_Mono(n_inputs,waitT,emissionT, emissionTR, epsT, cycles, interp_mode, 1);

% Free random couples
[Input,Response,tMax] = TI_Data_Multi_Mono(n_inputs,waitT,emissionT, emissionTR, epsT, cycles, interp_mode, 0);

% single images
[InputS,tMaxS] = TI_Data_Singular(n_inputs,waitT,emissionT, epsT, cycles, interp_mode);
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
%%
NuExtS = [InputS, Response]';
%%
NuExtTrans = NuExtAdj;
NuExtLearn = NuExtAdj;
NuExtAuto = NuExt;
%% Sets the stimulation protocol.

if OrthogonalStim
    Net.CParam.JExt = 2*(rand(Net.P,size(NuExtLearn, 1)) > 0.5) - 1; % Orthogonal (random spin) stimulation. (+-1)
    Net.CParam.JExt(:,1:n_inputs) = Net.CParam.JExt(:,1:n_inputs)*Options.StdJinput;
    Net.CParam.JExt(:,n_inputs+1:end) = Net.CParam.JExt(:,n_inputs+1:end)*Options.StdJback;
else
    Net.CParam.JExt = randn(Net.P,size(NuExtLearn, 1));
    Net.CParam.JExt(:,1:n_inputs) = Net.CParam.JExt(:,1:n_inputs)*Options.StdJinput;
    Net.CParam.JExt(:,n_inputs+1:end) = Net.CParam.JExt(:,n_inputs+1:end)*Options.StdJback;
end
%%
% add Directional Stimulous
Net.DirNoise.Intensity = NOISE_INT * Options.StdJinput;
Net.DirNoise.ENdx = n_inputs+1;
%% Numerically integrate the network dynamics
% Transient + Learning

[tIP, NuIP, Net] = RNNDynamicsDisc(Net, TransientPeriod, NuExtTrans, dt);
Net.t=0;
[tL, NuLP, Net] = RNNDynamicsDisc(Net, LearningPeriod, NuExtLearn, dt);

%%
% Plot inputs
alphabet = 'A':'Z';
colorsLearn = brewermap(n_inputs,'Paired');
colorsTrans = 0.5*colorsLearn;

%% Learn the coefficients "Wout"

OutOrig = NuExtLearn{1+n_inputs}(tL);

AddedNoise = NoiseForLearning*std(NuLP(:))*randn(size(NuLP));
X = NuLP+AddedNoise;
Wout = OutOrig/X;

OutL = Wout*NuLP;

%% Incorporates the learnt synaptic matrix and test the autonomous dynamics.

NetSelf = Net;
NetSelf.CParam.JOpen = NetSelf.CParam.J;
NetSelf.CParam.JToClose = Net.CParam.JExt(:,n_inputs+1:end)*Wout;
NetSelf.CParam.J = NetSelf.CParam.JOpen + NetSelf.CParam.JToClose;

NetSelf.MParam.nInputs = n_inputs;
NetSelf.MParam.waitT = waitT;
NetSelf.MParam.emissionT = emissionT;
NetSelf.MParam.emissionTR = emissionTR;
NetSelf.MParam.epsT = epsT;
NetSelf.MParam.cycles = cycles;
NetSelf.MParam.interp_mode = interp_mode;

NetSelf.CParam.JExt(:,n_inputs+1:end) = 0 * NetSelf.CParam.JExt(:,n_inputs+1:end);

AutonomousPeriod = AP*Options.Tau;
[tA, NuA, NetSelf] = RNNDynamicsDisc(NetSelf, AutonomousPeriod, NuExtAuto, dt);  
%%
% Read-out
OutA = Wout*NuA;
OutAOrig(1,:) = NuExtAuto{1+n_inputs}(tA);

%Plot reconstructed and expected output.
T = [tL tA];
if PlotFigures
    figure
    
    ax1 = subplot(2,3,1:3);
    hold on
    for k=1:n_inputs
        dist = 2.5*k;
        plot(tL/Options.Tau,dist+NuExtLearn{k}(tL),"Color",colorsLearn(k,:), "LineStyle","-",'LineWidth',1)
        plot(tA/Options.Tau,dist+NuExtAuto{k}(tA),"Color",colorsLearn(k,:), "LineStyle","-",'LineWidth',2)
        text(tA(end)/Options.Tau, dist, alphabet(k));
        plot(tL(end)/Options.Tau,dist+NuExtLearn{k}(tL(end)),"Marker","x", "Color",colorsLearn(k,:), "LineWidth",1.2, "MarkerSize",8)
    end
    
    text(0.01, -0.2, '$LearningPeriod$ $\rightarrow$',"Interpreter","latex",'Units','Normalized')
    text(tL(end)/tA(end), -0.2, '$|$ $AutonomousPeriod$ $\rightarrow$',"Interpreter","latex",'Units','Normalized')
    ylabel('Input')
    xlim([tL(1) tA(end)]/Options.Tau)
    ax1.XTick = [];
    ax1.YTick = [];
    ax1.XColor = [1 1 1];

    strTitle = ['N = ' num2str(Net.P),...
        '\hspace{0.2cm} g = ' num2str(Net.StdJinter),...
        '\hspace{0.2cm} $g_{in}$ = ' num2str(Options.StdJinput),...
        '\hspace{0.2cm} $||\beta||$ = ' num2str(norm(Wout), 3)];
    title(strTitle, "Interpreter","latex")

    %--
    ax2 = subplot(2,3,4:6);
    
    hold on
    plot(tL/Options.Tau,OutOrig(1,:),"Color",[0 0 0], "LineStyle","--",'LineWidth',1)
    plot(tA/Options.Tau,OutAOrig(1,:),"Color",[0 0 0], "LineStyle","--",'LineWidth',1)
    plot(tL/Options.Tau,OutL(1,:),"Color",[1 0 1], "LineStyle","-",'LineWidth',2)
    plot(tA/Options.Tau,OutA(1,:),"Color",[0 0 1], "LineStyle","-",'LineWidth',2)
    rectangle('Position',[tA(1)/Options.Tau -Options.th AutonomousPeriod/Options.Tau 2*Options.th], 'FaceColor',[0.5 0.5 0.5 0.5],'EdgeColor','w')
    
    xlim([tL(1) tA(end)]/Options.Tau)
    xlabel('Time, t/$\tau$', "Interpreter","latex")
    ylabel('Read-out output')        
    ax2.YTick = [];
    ylim([-1.5 1.5])

    linkaxes([ax1 ax2],'x');

end
%%
if PlotFigures
    a1 = 2200;
    a2 = 4100;
    figure
    ax1 = subplot(2,1,1);
    hold on
    for k=1:n_inputs
        dist = 2.6*k;
        plot(tA/Options.Tau,dist+NuExtAuto{k}(tA),"Color",colorsLearn(k,:), "LineStyle","-",'LineWidth',2)
        text(tA(end)/Options.Tau, dist, alphabet(k));
    end
    
    ylabel('Input')
    xlim([tA(a1) tA(a2)]/Options.Tau)
    ax1.XTick = [];
    ax1.YTick = [];
    ax1.XColor = [1 1 1];

    ax2 = subplot(2,1,2);
    
    hold on
    plot(tA/Options.Tau,OutA(1,:),"Color",[0 0 1], "LineStyle","-",'LineWidth',2)
    yticks(-10:0.5:10)
    xlim([tA(a1) tA(a2)]/Options.Tau)
    xlabel('Time, t/$\tau$', "Interpreter","latex")
    ylabel('Read-out output')        

    linkaxes([ax1 ax2],'x');

end
%%
a1 = 2200;
a2 = 4100;

%----
nun = 20;
colorsnu = brewermap(nun,'Set1');

figure
hold on
for k=1:nun
    plot(tA/Options.Tau, NuA(k,:)', color=colorsnu(k,:)*0.7)
end
xlim([tA(a1) tA(a2)]/Options.Tau)
ylabel('\nu')
xlabel('t/\tau')

%%
% select by symbolic distance
colorsSdist = brewermap(n_inputs-1,'Set3');

sdist = zeros(size(OutAOrig));
for k=1:n_inputs
    sdist = sdist + k*ceil(0.9*NuExtAuto{k} (tA)) - k*ceil(-0.9*NuExtAuto{k} (tA));
end
sdistNeg = sdist;
sdist = abs(sdist);

sdistMask = [];
for k=1:n_inputs-1
    sdistMask = [sdistMask {k}];
    sdistMask{k} = find(sdist == k);
end

%%
%=================================================================
% Learning
eventL = [];
counter = 0;
lastResp = 0;
tStart = tL(1);
minOutMax = inf;

relaxingCond = zeros(size(tL));
for t=tL
    if((NuExtLearn{end} (t) ~= 0) && (lastResp == 0)) % event starts
        tStart = t - DELAY;
        counter = counter + 1;
        eventL = [eventL {counter}];
        
        for k=1:n_inputs
            boolDx = ceil(0.9*NuExtLearn{k} (t));
            boolSx = ceil(-0.9*NuExtLearn{k} (t));
            
            if(boolDx)
                eventProp.image.dx = k;
            elseif(boolSx)
                eventProp.image.sx = k;
            end
        end
        
        
        if(t > waitT)
            relaxingCond = relaxingCond + ((tL < tStart) & (tL > tStart - 0.5*waitT));
        end
    elseif( ((NuExtLearn{end} (t) == 0) && (lastResp ~= 0)) || ... % event ends
            ((t == tL(end)) && (lastResp ~= 0)) )
        interval = find((tL >= tStart) & (tL <= t));
        eventProp.t = [tStart t];
        eventProp.Out = OutL(interval);
        [~,maxNdx] = max(abs(eventProp.Out));
        eventProp.maxOut = eventProp.Out(maxNdx);
        
        if(abs(eventProp.maxOut) < minOutMax)
            minOutMax = abs(eventProp.maxOut);
        end
        
        eventL{counter} = eventProp;
    end
    lastResp = NuExtLearn{end} (t);
    
end
%-----------------------------------------------------------
%% Event Definition
eventA = [];
counter = 0;
Sdist.counter = zeros(1,n_inputs-1); % num. events for each symb. distance
lastResp = 0;
eventProp.sdist = 0;
tStart = tA(1);
minOutMax = inf;

relaxingCond = zeros(size(tA));
for t=tA
    if((NuExtAuto{end} (t) ~= 0) && (lastResp == 0)) % eventA starts
        tStart = t - DELAY;
        counter = counter + 1;
        eventA = [eventA {counter}];
        
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
    elseif( ((NuExtAuto{end} (t) == 0) && (lastResp ~= 0)) || ... % eventA ends
            ((t == tA(end)) && (lastResp ~= 0)) )
        interval = find((tA >= tStart) & (tA <= t));
        eventProp.t = [tStart t];
        eventProp.Out = OutA(interval);
        [~,maxNdx] = max(abs(eventProp.Out));
        eventProp.maxOut = eventProp.Out(maxNdx);
        
        if(abs(eventProp.maxOut) < minOutMax)
            minOutMax = abs(eventProp.maxOut);
        end
        
        %5right_wrong = sign(eventProp.Out.*OutAOrig(interval) ); % element +1=right, -1=wrong 
        %acc = sum(right_wrong>0)/numel(interval);
        %eventProp.guess = (acc>0.5) ;
        
        eventA{counter} = eventProp;
    end
    lastResp = NuExtAuto{end} (t);
    
end
%%
%=======================================================================================
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
for k=1:numel(eventA)
    % maximum output
    Sdist.maxOut(eventA{k}.sdist) = Sdist.maxOut(eventA{k}.sdist) + abs(eventA{k}.maxOut);

    % reaction times
    interval = find((tA >= eventA{k}.t(1)) & (tA <= eventA{k}.t(end)));
    tInt = tA(interval);
    OutInt = OutA(interval);
    OutOrigInt = OutAOrig(interval);
    
    guess = 0;
    if(~isempty(find(abs(OutInt) > th)))
        ndxOTH = find(abs(OutInt) > th);
        tTH = tInt(ndxOTH);
        Sdist.reactTimes(eventA{k}.sdist) = Sdist.reactTimes(eventA{k}.sdist) + (tTH(1)-eventA{k}.t(1)); 
    
        %accuracy
        guess = (sign(OutInt(ndxOTH(1))*OutOrigInt(ndxOTH(1))) > 0);

        counter_oTH_Heat(eventA{k}.image.dx,eventA{k}.image.sx) = counter_oTH_Heat(eventA{k}.image.dx,eventA{k}.image.sx) + 1;
        SD_oTH_counter(eventA{k}.sdist) = SD_oTH_counter(eventA{k}.sdist) + 1; %just for reaction times
    end
    Sdist.acc(eventA{k}.sdist) = Sdist.acc(eventA{k}.sdist) + guess;

    % accuracy heatmap
    AccHeat(eventA{k}.image.dx,eventA{k}.image.sx) = AccHeat(eventA{k}.image.dx,eventA{k}.image.sx) + guess;
    counterHeat(eventA{k}.image.dx,eventA{k}.image.sx) = counterHeat(eventA{k}.image.dx,eventA{k}.image.sx) + 1;

end

SymDist.Acc = Sdist.acc;
SymDist.MaxOut = Sdist.maxOut ./ Sdist.counter;
SymDist.ReactTimes = Sdist.reactTimes ./ SD_oTH_counter;

if(EFFECT_ONLY_OTH)
    SymDist.Acc = SymDist.Acc ./ SD_oTH_counter;
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
%=======================================================================================
%%
colorNoise = [140 240 220]/255;
% SDE EFFECT
if PlotFigures
    figure
    hold on
    %title('SDE ACC')

    f = griddedInterpolant(1:n_inputs-1,SymDist.Acc,'pchip');
    x = 1:0.1:(n_inputs-1);
    plot(x,f(x),"Color",colorNoise,"LineWidth", 3,"HandleVisibility","off")
    plot(1:n_inputs-1, SymDist.Acc,"LineStyle","none", "LineWidth", 4,...
        "Marker", "o", "MarkerSize", 12, "Color", colorNoise,"MarkerFaceColor","w")   

    ylim([0.0 1])
    ylabel('Accuracy')
    xlabel('Symbolic Distance')
end
%%
if PlotFigures
    figure
    hold on

    f = griddedInterpolant(1:n_inputs-1,SymDist.ReactTimes,'pchip');
    x = 1:0.1:(n_inputs-1);
    plot(x,f(x),"Color",colorNoise,"LineWidth", 3,"HandleVisibility","off")
    plot(1:n_inputs-1, SymDist.ReactTimes,"LineStyle","none", "LineWidth", 4,...
        "Marker", "o", "MarkerSize", 12, "Color", colorNoise,"MarkerFaceColor","w")   


    ylabel('reactTimes')
    xlabel('Symbolic Distance')
    %title('SDE RT')
end
%%
if PlotFigures
    figure
    plot(AnchorAcc,"LineStyle","-","Marker","o","MarkerFaceColor",'w',"MarkerSize",8,...
        "MarkerEdgeColor",[0.1 0.5 0.3],'LineWidth',2, 'Color', [0.1 0.5 0.3])
    ylim([0.0,1])
    xticks(1:n_inputs)
    xticklabels(alphabet(1:n_inputs)')
    xlim([1,(n_inputs)])
    xlabel('Symbol')
    ylabel('Mean Accuracy')
end
%%
if PlotFigures
    figure
    plot(1:n_inputs-1, SymDist.MaxOut,"LineStyle","-","Marker","o","MarkerFaceColor",'w',"MarkerSize",8,...
        "MarkerEdgeColor",[0 0 0],"Color",[0 0 0],'LineWidth',1.5);
    ylabel('maxOutput')
    xlabel('Symbolic Distance')
    title('Max Output vs SDist')
end
%%
corrTh = 0.4;
neuroPercs = zeros(1,n_inputs-1);
neuroMean = zeros(1,n_inputs-1);
for k=1:n_inputs-1
    M = [OutA(sdist==k);NuA(:,sdist==k)]';
    Rm = corrcoef(M);
    Rc = Rm(1,2:end);
    indxM_OvTh = (abs(Rc)>corrTh);
    neuroPercs(k) = sum(indxM_OvTh)/Net.P;
    neuroMean(k) = mean(abs(Rc));
end
%%
if PlotFigures
    figure
    plot(neuroPercs,"LineStyle","-","Marker","o","MarkerFaceColor",'w',"MarkerSize",8,...
        "MarkerEdgeColor",[0.9 0.8 0.1],'LineWidth',1.5, 'Color', [0.9 0.8 0.1])
    hold on
    plot(neuroMean)
    ylim([0.0,1])
    xlim([1,(n_inputs-1)])
    xlabel('Symb. Distance')
    ylabel('CorOvTH (%)')
end
%%
% Accuracy Grid
if PlotFigures
    figure
    heatmap(alphabet(1:n_inputs)',alphabet(1:n_inputs)',AccHeat);
    title('Accuracy Table E')
end

%%
signIndxL_R = find(OutOrig>0); % only right solution
signIndxA_R = find(OutA>0);

signIndxL_L = find(OutOrig<0); % only left solution
signIndxA_L = find(OutA<0);
%---------------------------------
[PCsA,~,EVs_A] = pca(NuA');
PrjsL_A = NuLP' * PCsA;
PrjsA_A = NuA' * PCsA;

%%
if(PlotFigures)
    figure
    eVar = cumsum(EVs_A)/sum(EVs_A)*100;
    plot(eVar,'ko-','MarkerFaceColor','w');
    set(gca,'TickDir','out','YLim',[0 100],'XLim',[1 Net.P],'XScale','log');
    xlabel('PCs')
    ylabel('Explained variance (%)')
    title('Autonomous')

    figure
    
    subplot(1,2,1) 
    plot(PrjsL_A(:,1),PrjsL_A(:,2),'--', 'color', 0.8*[1 1 1])
    hold('on')
    
    scatter(PrjsL_A(:,1),PrjsL_A(:,2), 'MarkerEdgeColor',[0 0 0], "MarkerEdgeAlpha",0.2)
    
    %---
    colorSign = [
        0.8 0 0.5;
        0.2 0.5 0.7
        ];
    plot(PrjsL_A(signIndxL_R,1),PrjsL_A(signIndxL_R,2),'.', 'color', colorSign(1,:))
    plot(PrjsL_A(signIndxL_L,1),PrjsL_A(signIndxL_L,2),'.', 'color', colorSign(2,:))
    
    colormap(gca,colorSign);
    cb = colorbar('north', 'Ticks', [0.25 0.75],...
        'TickLabels', ['Right'; ' Left'], 'FontSize',8);
    cb.Label.FontSize = 6;
    ax = gca;
    axpos = ax.Position; 
    cb.Position(4) = 0.5*cb.Position(4);
    cb.Position(2) = 1.05*cb.Position(2);
    ax.Position = axpos;
    
    title('Learning');
    set(gca,'TickDir','out');
    xlabel('PC_1')
    ylabel('PC_2')
    
    %---------------------------------
    
    subplot(1,2,2)
    plot(PrjsA_A(:,1),PrjsA_A(:,2),'--', 'Color',0.8*[1 1 1])
    hold on
    scatter(PrjsA_A(:,1),PrjsA_A(:,2), 'MarkerEdgeColor',[0 0 0], "MarkerEdgeAlpha",0.2)
    
    
    plot(PrjsA_A(signIndxA_R,1),PrjsA_A(signIndxA_R,2),'.', 'color', colorSign(1,:))
    plot(PrjsA_A(signIndxA_L,1),PrjsA_A(signIndxA_L,2),'.', 'color', colorSign(2,:))
    
    colormap(gca,colorSign);
    cb = colorbar('north', 'Ticks', [0.25 0.75],...
        'TickLabels', [' Left'; 'Right'], 'FontSize',8);
    cb.Label.FontSize = 6;
    ax = gca;
    axpos = ax.Position; 
    cb.Position(4) = 0.5*cb.Position(4);
    cb.Position(2) = 1.05*cb.Position(2);
    ax.Position = axpos;
    
    title('Autonomous');
    set(gca,'TickDir','out');
    xlabel('PC_1')
end
%%
imgIndxL = [];
imgIndxA = [];
imgDxIndxA = cell(1,N_INPUTS);
for k=1:n_inputs
    imgDxIndxATemp = find((NuExtAuto{k} (tA)) > 0);
    imgDxIndxA{k} = imgDxIndxATemp(1:2:end);
end

if PlotFigures
    figure
    plot3(PrjsA_A(:,1),PrjsA_A(:,2), PrjsA_A(:,3), 'Color', 0.7*[1 1 1], "LineStyle","--")
    hold on
    for k=1:n_inputs
        plot3(PrjsA_A(imgDxIndxA{k},1),PrjsA_A(imgDxIndxA{k},2), PrjsA_A(imgDxIndxA{k},3),...
            '.','Color', colorsLearn(k,:), 'MarkerSize',40, "LineWidth",2)
        hold on
    end
    grid on
    colormap(colorsLearn);
    cb = colorbar('northoutside', 'TickLabels', {});
    ax = gca;
    axpos = ax.Position; 
    cb.Position(4) = 0.5*cb.Position(4);
    cb.Position(2) = cb.Position(2)+0.085;
    ax.Position = axpos;
end
%%
colorSignMulti = flip([
    206 0 127  ;
    190 54 135 ;
    176 73 142 ;
    162 86 148 ;
    149 95 154 ;
    137 102 159;
    125 108 163;
    113 113 166;
    102 116 169;
    90 120 172 ;
    74 123 175 ;
    50 127 179  
    ]/255);

xlims = [min(PrjsA_A(:,1))-0.5 max(PrjsA_A(:,1))+0.5];
if PlotFigures
    figure
    subplot(8,1,1:3)
    plot(PrjsA_A(:,1),PrjsA_A(:,2),"LineStyle","none", ...
        "Marker","o", "MarkerEdgeColor",0.88*[1 1 1], "MarkerSize", 5)
    hold on
    plot(PrjsA_A(:,1),PrjsA_A(:,2),"LineStyle","none", ...
        "Marker",".", "MarkerEdgeColor",0.3*[1 1 1], "MarkerSize", 3)
    
    
    for k=1:n_inputs-1
    scatter(PrjsA_A(sdistMask{k},1),PrjsA_A(sdistMask{k},2),...
           50,'MarkerEdgeColor',0.3*[1 1 1],'Marker',".","MarkerEdgeAlpha",0.3)
    end
    for k=n_inputs-1:-1:1
        ndxU = (OutA(sdistMask{k}) > 0);
        ndxD = (OutA(sdistMask{k}) < 0);
        A = PrjsA_A(sdistMask{k},1);
        B = PrjsA_A(sdistMask{k},2);
        scatter(A(ndxU),B(ndxU),...
           50,'MarkerEdgeColor',colorSignMulti(n_inputs-1+k,:),...
        "MarkerEdgeAlpha",0.25,'Marker',"o")
        scatter(A(ndxD),B(ndxD),...
           50,'MarkerEdgeColor',colorSignMulti(n_inputs-k,:),...
        "MarkerEdgeAlpha",0.25,'Marker',"o")    
    end
        
    
    ylabel("PC_2")
    set(gca,'TickDir','out','XColor','w');
    box off
    h = gca;
    h.XAxis.Visible = 'off';
    xlim(xlims)
    
    
    %------------------
    subplot(8,1,4:8)
    hold on

    plot(PrjsA_A(:,1),OutA,"LineStyle","none", ...
        "Marker","o", "MarkerEdgeColor",0.92*[1 1 1], "MarkerSize", 5)
    hold on
    scatter(PrjsA_A(:,1),OutA,2,...
        "Marker",".", "MarkerEdgeColor",0.3*[1 1 1])

    
    for k=1:n_inputs-1
    scatter(PrjsA_A(sdistMask{k},1),OutA(sdistMask{k}),...
           60,'MarkerEdgeColor',0.3*[1 1 1],'Marker',".")
    hold on
    end
    
    for k=n_inputs-1:-1:1
        ndxU = (OutA(sdistMask{k}) > 0);
        ndxD = (OutA(sdistMask{k}) < 0);
        A = PrjsA_A(sdistMask{k},1);
        B = OutA(sdistMask{k});
        scatter(A(ndxU),B(ndxU),...
           100,'MarkerEdgeColor',colorSignMulti(n_inputs-1+k,:),...
        "MarkerEdgeAlpha",0.1,'Marker',"o")
        scatter(A(ndxD),B(ndxD),...
           100,'MarkerEdgeColor',colorSignMulti(n_inputs-k,:),...
        "MarkerEdgeAlpha",0.1,'Marker',"o")    
    end
    
    
    xlim(xlims)
    xlabel('PC_1')

    yticks(-N_INPUTS:N_INPUTS)
    set(gca,'Layer','top','TickDir','out','Box','off')
    ylabel('Output')
    
end
%%


colorSignMulti = flip([
    206 0 127  ;
    190 54 135 ;
    176 73 142 ;
    162 86 148 ;
    149 95 154 ;
    137 102 159;
    125 108 163;
    113 113 166;
    102 116 169;
    90 120 172 ;
    74 123 175 ;
    50 127 179  
    ]/255);

if PlotFigures
    figure
    
    plot(PrjsA_A(:,1),OutA,"LineStyle","none", ...
        "Marker",".", "MarkerEdgeColor",0.92*[1 1 1])
    
    hold on

    for k=1:n_inputs-1
    scatter(PrjsA_A(sdistMask{k},1),OutA(sdistMask{k}),...
           100,'MarkerEdgeColor',0.3*[1 1 1],'Marker',"o")
    end
    
    for k=n_inputs-1:-1:1
        ndxU = (OutA(sdistMask{k}) > 0);
        ndxD = (OutA(sdistMask{k}) < 0);
        A = PrjsA_A(sdistMask{k},1);
        B = OutA(sdistMask{k});
        scatter(A(ndxU),B(ndxU),...
           150,'MarkerEdgeColor',colorSignMulti(n_inputs-1+k,:),...
        "MarkerEdgeAlpha",0.25,'Marker',".")
        scatter(A(ndxD),B(ndxD),...
           150,'MarkerEdgeColor',colorSignMulti(n_inputs-k,:),...
        "MarkerEdgeAlpha",0.25,'Marker',".")    
    end
    
    
    
    colormap(colorSignMulti);
    cbTicks = (0.5:2*n_inputs)/(2*(n_inputs-1));
    ss = zeros(1,2*(n_inputs-1));
    ss(1:(n_inputs-1)) = -n_inputs+1:-1;
    ss(n_inputs:end) = 1:n_inputs-1;
    cb = colorbar('eastoutside', 'Ticks', cbTicks,...
        'TickLabels', ss, 'FontSize',8,"TickDirection","out");
    cb.Label.String = 'Symb. Dist.';
    cb.Label.FontSize = 8;
    ax = gca;
    axpos = ax.Position; 
    cb.Position(3) = 0.5*cb.Position(3);
    ax.Position = axpos;
    
    yticks(-N_INPUTS:N_INPUTS)
    set(gca,'Layer','top','TickDir','out','Box','off')
    xlabel('PC_1')
    ylabel('Read-out Output')
end
%%
eW = Wout/norm(Wout);
exVarW = var(eW * NuA) / sum(EVs_A)
%%
if(PlotFigures)
    figure
    
    plot3(PrjsA_A(:,1),PrjsA_A(:,2),PrjsA_A(:,3),'-', 'Color', 0.7*[1 1 1],"LineWidth",1.5)
    
    hold on
    
    plot3(PrjsA_A(signIndxA_R,1),PrjsA_A(signIndxA_R,2),PrjsA_A(signIndxA_R,3),'.', ...
        'Color', colorSign(1,:), 'MarkerSize',15)
    plot3(PrjsA_A(signIndxA_L,1),PrjsA_A(signIndxA_L,2),PrjsA_A(signIndxA_L,3),'.', ...
        'Color', colorSign(2,:), 'MarkerSize', 15)
    colormap(gca,colorSign);
    cb2 = colorbar('northoutside', 'Ticks', [0.25 0.75], ...
        'TickLabels', ['Right'; ' Left'],'FontSize',10,"TickDirection","out");
    ax2 = gca;
    axpos = ax2.Position; 
    cb2.Position(4) = 0.5*cb2.Position(4);
    cb2.Position(3) = 0.5*cb2.Position(3);
    cb2.Position(2) = cb2.Position(2)-0.15;
    cb2.Position(1) = cb2.Position(1)+0.2;
    ax2.Position = axpos;
    grid on
    xlabel('PC_1')
    ylabel('PC_2')
    zlabel('PC_3')
end
%% Linear Decoder

% Learning
MaskCoupleL = zeros(n_inputs, numel(tL));
count = 1;
for t = 1:numel(tL)
    if((tL(t) > eventL{count}.t(1)) && (tL(t) < eventL{count}.t(2)))
        MaskCoupleL(eventL{count}.image.dx,t) = +1;
        MaskCoupleL(eventL{count}.image.sx,t) = -1;
    elseif((tL(t) == eventL{count}.t(2)) && (count < numel(eventL)))
        count = count +1;
    end
end
%------------------------------------------
% Autonomous
MaskCoupleA = zeros(n_inputs, numel(tA));
count = 1;
for t = 1:numel(tA)
    if((tA(t) > eventA{count}.t(1)) && (tA(t) < eventA{count}.t(2)))
        MaskCoupleA(eventA{count}.image.dx,t) = 1;
        MaskCoupleA(eventA{count}.image.sx,t) = -1;
    elseif((tA(t) == eventA{count}.t(2)) && (count < numel(eventA)))
        count = count +1;
    end
end

S_L = MaskCoupleL*pinv(X);

%% Singular Stimulation

SingPeriod = SP*Options.Tau;
[tS, NuS, ~] = RNNDynamicsDisc(NetSelf, SingPeriod, NuExtS, dt);  

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
        interval = find((tS >= tStart) & (tS <= t));
        eventProp.t = [tStart t];
        eventProp.Out = OutRank(interval);
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
if PlotFigures
    figure
    
    ax1 = subplot(2,3,1:3);
    hold on
    for k=1:n_inputs
        dist = 2.5*k;
        plot(tS/Options.Tau,dist+NuExtS{k}(tS),"Color",colorsLearn(k,:), "LineStyle","-",'LineWidth',2)
        text(tS(end)/Options.Tau, dist, alphabet(k));
    end
    
    ylabel('Input')
    xlim([tS(1) tS(end)]/Options.Tau)
    ax1.XTick = [];
    ax1.YTick = [];
    ax1.XColor = [1 1 1];


    %--
    ax2 = subplot(2,3,4:6);
    plot(tS/Options.Tau,OutRank(1,:),"Color",[0.2 0.4 1], "LineStyle","-",'LineWidth',2)
    hold on
    for k=1:numel(event)
        plot(event{k}.t/Options.Tau,event{k}.maxOut,"LineStyle","none","Marker",".","MarkerFaceColor",[1 0 0],"MarkerSize",6,...
            "MarkerEdgeColor",[1 0 0.2]);
    end
    xlim([tS(1) tS(end)]/Options.Tau)
    xlabel('Time, t/$\tau$', "Interpreter","latex")
    ylabel('Read-out output')        
    ax2.YTick = [];
        
    linkaxes([ax1 ax2],'x');
end
%%
distrImg = nan(n_inputs, numel(event));
for k=1:numel(event)
    distrImg(event{k}.image.dx,k) = event{k}.maxOut;
end
%%
YL = 1.5;
if PlotFigures
    figure
    for k=1:n_inputs
        histogram(distrImg(k,~isnan(distrImg(k,:))),30,"Normalization","pdf","FaceColor",colorsLearn(k,:),"FaceAlpha",0.1,"EdgeColor",colorsLearn(k,:),"EdgeAlpha",0.2)
        hold on
        pd = fitdist(distrImg(k,:)','Normal');
        xgrid = linspace(min(distrImg(k,:)),max(distrImg(k,:)),100)';
        pdfEst = pdf(pd,xgrid); 
        line(xgrid,pdfEst,'Color',colorsLearn(k,:),'LineWidth',2)
        plot(mean(distrImg(k,~isnan(distrImg(k,:)))), 0, 'LineStyle', 'none', 'Marker', '|', ...
            'Color',colorsLearn(k,:),'MarkerSize',10,'MarkerFaceColor',colorsLearn(k,:),'LineWidth',4)
    end
    colormap(gca,colorsLearn);
    xlabel('Mental Projection')
    cbTicks = (0.5:n_inputs+1)/(n_inputs);
    cb = colorbar('north', 'Ticks', cbTicks,...
        'TickLabels', alphabet(1:n_inputs)', 'FontSize',12);
    ax = gca;
    axpos = ax.Position; 
    cb.Position(4) = 0.5*cb.Position(4);
    cb.Position(2) = 1.05*cb.Position(2);
    ax.Position = axpos;
    ax.YColor = [1 1 1];
    ax.YTick = [];
    box off
    ylim([0 YL]);

end
%%
if PlotFigures
    figure
    for k=1:n_inputs
        histogram(distrImg(k,~isnan(distrImg(k,:))),30,"Normalization","pdf","FaceColor",colorsLearn(k,:),"FaceAlpha",0.1,"EdgeColor",colorsLearn(k,:),"EdgeAlpha",0.2)
        hold on
        pd = fitdist(distrImg(k,:)','kernel');
        xgrid = linspace(min(distrImg(k,:)),max(distrImg(k,:)),100)';
        pdfEst = pdf(pd,xgrid); 
        line(xgrid,pdfEst,'Color',colorsLearn(k,:),'LineWidth',2)
        plot(mean(distrImg(k,~isnan(distrImg(k,:)))), 0, 'LineStyle', 'none', 'Marker', '|', ...
            'Color',colorsLearn(k,:),'MarkerSize',10,'MarkerFaceColor',colorsLearn(k,:),'LineWidth',4)
    end
    colormap(gca,colorsLearn);
    xlabel('Mental Projection')
    cbTicks = (0.5:n_inputs+1)/(n_inputs);
    cb = colorbar('north', 'Ticks', cbTicks,...
        'TickLabels', alphabet(1:n_inputs)', 'FontSize',12);
    %cb.Label.String = 'Symb.';
    %cb.Label.FontSize = 10;
    ax = gca;
    axpos = ax.Position; 
    cb.Position(4) = 0.5*cb.Position(4);
    cb.Position(2) = 1.05*cb.Position(2);
    ax.Position = axpos;
    ax.YColor = [1 1 1];
    ax.YTick = [];
    %title("Max Out Points - Event Counts per Image: " + num2str(C))
    %grid on
    box off
    ylim([0 YL]);

end
%%
RM = round(0.5*(numel(event{1}.Out) - round(0.1*numel(event{1}.Out)))); %select the central 10% of the event
distrImgContinous = nan(n_inputs, numel(event)*numel(event{1}.Out));
for k=1:numel(event)
    centralSignal = event{k}.Out(RM:end-RM);
    distrImgContinous(event{k}.image.dx,((k-1)*numel(centralSignal)+1):k*numel(centralSignal)) = centralSignal;
end
%%
if PlotFigures
    figure
    for k=1:n_inputs
        histogram(distrImgContinous(k,~isnan(distrImgContinous(k,:))),100,"Normalization","pdf", ...
            "FaceColor",colorsLearn(k,:),"FaceAlpha",0.1,"EdgeColor",colorsLearn(k,:),"EdgeAlpha",0.2);
        hold on
        pd = fitdist(distrImgContinous(k,:)','Normal');
        xgrid = linspace(min(distrImgContinous(k,:)),max(distrImgContinous(k,:)),100)';
        pdfEst = pdf(pd,xgrid); 
        line(xgrid,pdfEst,'Color',colorsLearn(k,:),'LineWidth',2);
        plot(mean(distrImgContinous(k,~isnan(distrImgContinous(k,:)))), 0, 'LineStyle', 'none', 'Marker', '|', ...
            'Color',colorsLearn(k,:),'MarkerSize',10,'MarkerFaceColor',colorsLearn(k,:),'LineWidth',4);
    end
    colormap(gca,colorsLearn);
    xlabel('Mental Projection')
    cbTicks = (0.5:n_inputs+1)/(n_inputs);
    cb = colorbar('north', 'Ticks', cbTicks,...
        'TickLabels', alphabet(1:n_inputs)', 'FontSize',12);
    ax = gca;
    axpos = ax.Position; 
    cb.Position(4) = 0.5*cb.Position(4);
    cb.Position(2) = 1.05*cb.Position(2);
    ax.Position = axpos;
    ax.YColor = [1 1 1];
    ax.YTick = [];
    box off
    ylim([0 YL]);

end
%%
if PlotFigures
    figure
    for k=1:n_inputs
        histogram(distrImgContinous(k,~isnan(distrImgContinous(k,:))),100,"Normalization","pdf", ...
            "FaceColor",colorsLearn(k,:),"FaceAlpha",0.1,"EdgeColor",colorsLearn(k,:),"EdgeAlpha",0.2)
        hold on
        pd = fitdist(distrImgContinous(k,:)','kernel');
        xgrid = linspace(min(distrImgContinous(k,:)),max(distrImgContinous(k,:)),100)';
        pdfEst = pdf(pd,xgrid); 
        line(xgrid,pdfEst,'Color',colorsLearn(k,:),'LineWidth',2);
        plot(mean(distrImgContinous(k,~isnan(distrImgContinous(k,:)))), 0, 'LineStyle', 'none', 'Marker', '|', ...
            'Color',colorsLearn(k,:),'MarkerSize',10,'MarkerFaceColor',colorsLearn(k,:),'LineWidth',4)
    end
    colormap(gca,colorsLearn);
    xlabel('Mental Projection')
    box off
    cbTicks = (0.5:n_inputs+1)/(n_inputs);
    cb = colorbar('north', 'Ticks', cbTicks,...
        'TickLabels', alphabet(1:n_inputs)', 'FontSize',12);
    ax = gca;
    axpos = ax.Position; 
    cb.Position(4) = 0.5*cb.Position(4);
    cb.Position(2) = 1.05*cb.Position(2);
    ax.Position = axpos;
    ax.YColor = [1 1 1];
    ax.YTick = [];
    ylim([0 YL]);

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
if PlotFigures
    figure
    heatmap(alphabet(1:n_inputs)',alphabet(1:n_inputs)',wrs_pv)
    title('Wilcoxon Test on Max Distr.')
end
%%
if PlotFigures
    figure
    plot(wrs_im)
    xlabel('symbol')
    ylabel('mean wilcoxon pv')
end
%%
if PlotFigures
    figure
    plot(wrsC_im)
    xlabel('symbol')
    ylabel('mean wilcoxonC pv')
end
%%
toc
