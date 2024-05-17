%% 
clear all
close all
%% Create mean vectors (Symbols)
M = 7; % #symbols
N = 100; % dimension
S = randn(N,M);

alp = ('A':'Z')';
%voc = alp(1:M);
voc = alp(M:-1:1);
SAVE_FIGURES = false;

%normalize vectors
for k=1:M
    S(:,k) = S(:,k)/norm(S(:,k));
end

%% Orthogonalizing
[ST,~] = mgson(S);

%% Adding Noise on couples
N_TrialsL = 10000; % trials per couple
NoiseIntensity = 0.0;
sigma = 0.0;

Xnoise = nan(N,N_TrialsL*(2*(M-1)));
Onoise = ones(1,size(Xnoise,2));
c = 1;
for k=1:(M-1)
    for v=1:N_TrialsL
        Snoise = ST + NoiseIntensity*randn(N,M);
        vec = Snoise(:,k+1)-Snoise(:,k);
        Xnoise(:,c) = vec;
        Onoise(c) = sign(log(1+1/k)+sigma*randn());
        c = c+1;
    end
end
Xnoise(:,N_TrialsL*(M-1)+1:end) = -Xnoise(:,1:N_TrialsL*(M-1));

Onoise(N_TrialsL*(M-1)+1:end) = -Onoise(1:N_TrialsL*(M-1));

XiNoisy_ti = (Onoise*pinv(Xnoise))';

%% violin plot

sigmaS = 0.3/sqrt(2);
Y = zeros(N_TrialsL, M);
for n=1:N_TrialsL
    Y(n,:) = (ST+sigmaS*randn(N,M))'*XiNoisy_ti;
end
mY = mean(Y);

%norm(XiNoisy_ti')
%norm(XiNoisy_ti' - sum((ST'*XiNoisy_ti) .* ST'))
%-----------
colors = [255 51 51;
    250 201 56;
    168 245 61;
    45 211 86;
    139 243 243;
    76 114 230;
    153 81 225]/255;
colors = flipud(colors);

alph = 0.3;
YLIM = [-10, +8];
figure
%scatter(OSnoise_T, all_proj, 150,'k.')
violin(Y,'facecolor',colors,'edgecolor','k',...
'bw',0.3,...
'mc','k',...
'linewidth', 5,...
'plotlegend', false,...
'medc','',14)

hold on
for k=1:M
    scatter(k, mY(1)+(k-1), 2000,'w.')
    scatter(k, mY(1)+(k-1), 1000,'.', 'MarkerEdgeColor',colors(k,:))
end
xticks(1:M)
xticklabels(voc)
xlim([0.5,M+0.5])
ylim([YLIM(1),YLIM(2)])
ylabel("projection")
xlabel('Items')
title("sigma = "+num2str(sigma))
%xticklabels(voc)
ax = gca;
set(ax, 'TickDir', 'out', 'Box', 'off')
pbaspect([1.0 1.5 1])
if SAVE_FIGURES
    print('-painters','-dpdf',"Rew_proj_"+num2str(sigma) + ".pdf")
end

%% Test on all couples
N_TrialsD = 200;
%NoiseIntensity = 0.3;
NoiseIntensity = 0.0;
XnoiseSD = nan(N,M*(M-1)*N_TrialsD);
OnoiseSD = nan(1,M*(M-1)*N_TrialsD);
ImageNoiseR = nan(1,M*(M-1)*N_TrialsD);
ImageNoiseL = nan(1,M*(M-1)*N_TrialsD);
counter = 1;
for d=1:M
    for s=1:M
        sd = (s-d); %sdist
        if(sd~=0)
            for v=1:N_TrialsD
                Snoise = ST + NoiseIntensity*randn(N,M);
                vec = Snoise(:,s)-Snoise(:,d);
                XnoiseSD(:,counter) = vec;
                OnoiseSD(counter) = sd;
                ImageNoiseR(counter) = d;
                ImageNoiseL(counter) = s;
                counter = counter+1;
            end
        end
    end
end
%%
Outs = XnoiseSD'*XiNoisy_ti;
figure
yline(0, '-', LineWidth=2, Color=0.7*[1 1 1], Alpha=0.5)
hold on
plot([-(M+0.5) (M+0.5)],[-(M+0.5) (M+0.5)],'--r', LineWidth=2)
scatter(OnoiseSD,Outs, 100, 'k.')
xlim([-(M+0.5) (M+0.5)])
xticks([[-M+1:-1],[1:M-1]])
xlabel("SD")
ylabel("Readout")
title("Single Param")
ax = gca;
set(ax, 'TickDir', 'out')
if SAVE_FIGURES
    hgexport(gcf,'XiProb_TI_proj');
end

%% Accuracy check
RespXi = sign(XiNoisy_ti'*XnoiseSD);
Sd_mask = abs(OnoiseSD);
Resp = sign(OnoiseSD);

Acc = zeros(1,M-1);
for k=1:(M-1)
    ndx = find(Sd_mask == k);
    Acc(k) = sum(Resp(ndx)==RespXi(ndx)) / numel(ndx);
end

AccS = zeros(1,M);
for k=1:M
    ndx = find((ImageNoiseR == k) | (ImageNoiseL == k));
    AccS(k) = sum(Resp(ndx)==RespXi(ndx)) / numel(ndx);
end
%%
figure
plot(Acc, 'o-k', LineWidth=2.5, MarkerSize=8, MarkerFaceColor='white')
ylim([0.5 1])
xlim([1,M-1])
ylabel("Accuracy")
xlabel("|SD|")
title("Single Param")
ax = gca;
set(ax, 'TickDir', 'out')
pbaspect([1 1.5 1])
if SAVE_FIGURES
    hgexport(gcf,'XiProb_TI_Acc');
end
%%
figure
plot(AccS, 'o-k', LineWidth=3, MarkerSize=8, MarkerFaceColor='white')
ylim([0.5 1])
ylabel("Accuracy")
title("Single Param")
xticklabels(voc)
xlim([1,M])
ax = gca;
set(ax, 'TickDir', 'out')
pbaspect([1 1.5 1])
if SAVE_FIGURES
    hgexport(gcf,'XiProb_TI_Anch');
end

%% Testing different parameters
N = 100;
M=7;
sigmas = [0.0, 0.25, 0.5];
N_TrialsL=10000;
N_TrialsT=10000;
noiseL=0.0;
noiseT=0.3;

projs = zeros(M,numel(sigmas));
Accs=zeros(M-1,numel(sigmas));
AccSs=zeros(M,numel(sigmas));
for k=1:numel(sigmas)
    [projs(:,k),Accs(:,k),AccSs(:,k)] = GML_ProbRewards(N,M,sigmas(k),N_TrialsL,N_TrialsT,noiseL,noiseT);
end
%%
SAVE_FIGURES = false;
alp = ('A':'Z')';
voc = alp(1:M);
hue = flip(brewermap(numel(sigmas)+3, 'Greys'));
hue = flip(hue(1:numel(sigmas)+0,:));

figure
hold on
for k=1:numel(sigmas)
    plot(projs(:,k), 'o-k', LineWidth=2.5, MarkerSize=8, MarkerFaceColor='white', color=hue(k,:))
end
xticklabels(voc)
xlim([1 M])
grid off
ylabel('Projection')
ax = gca;
set(ax, 'TickDir', 'out')
pbaspect([1 1.5 1])
if SAVE_FIGURES
    hgexport(gcf,'XiProb_proj');
end
%%
figure
hold on
for k=1:numel(sigmas)
    plot(Accs(:,k), 'o-k', LineWidth=2.5, MarkerSize=8, MarkerFaceColor='white', color=hue(k,:))
end
xlim([1 M-1])
ylim([0.5 1.0])
grid off
ylabel('accuracy')
ax = gca;
set(ax, 'TickDir', 'out')
pbaspect([1 1.5 1])
if SAVE_FIGURES
    hgexport(gcf,'XiProb_TI_Acc');
end
%%
figure
hold on
for k=1:numel(sigmas)
    plot(AccSs(:,k), 'o-k', LineWidth=2.5, MarkerSize=8, MarkerFaceColor='white', color=hue(k,:))
end
xticklabels(voc)
xlim([1 M])
ylim([0.5,1.0])
grid off
ylabel('mean accuracy')
ax = gca;
set(ax, 'TickDir', 'out')
pbaspect([1 1.5 1])
if SAVE_FIGURES
    hgexport(gcf,'XiProb_TI_Anch');
end
