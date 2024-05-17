%% 
clear all
close all
%% Create mean vectors (Symbols)
M = 7; % #symbols
N = 100; % dimension
S = randn(N,M);

alp = ('A':'Z')';
voc = alp(1:M);
SAVE_FIGURES = false;

%normalize vectors
for k=1:M
    S(:,k) = S(:,k)/norm(S(:,k));
end
%% Orthogonalizing
[ST,~] = mgson(S);

%% Adding different num. of trials on single item
NNET = 100;
projs = zeros(NNET,M);
mt = 5;
N_TrialsLs = 50:-mt:(50-(M-1)*mt);
Bias = 0;

% crea funzione in 
NoiseIntensity = 0.3;

for n=1:NNET
    Xnoise_UN = nan(N,sum(N_TrialsLs));
    Onoise_UN = ones(1,size(Xnoise_UN,2));
    c = 1;
    for k=1:M
        for v=1:N_TrialsLs(k)
            Snoise = ST + NoiseIntensity*randn(N,M);
            Xnoise_UN(:,c) = Snoise(:,k);
            Onoise_UN(c) = k+Bias;
            c = c+1;
        end
    end
    
    XiNoisy_UN = (Onoise_UN*pinv(Xnoise_UN))';
    projs(n,:) = ST'*XiNoisy_UN;
end
%%
mt2 = 2;
projs2 = zeros(NNET,M);
N_TrialsLs2 = 50:-mt2:(50-(M-1)*mt2)
for n=1:NNET
    Xnoise_UN = nan(N,sum(N_TrialsLs2));
    Onoise_UN = ones(1,size(Xnoise_UN,2));
    c = 1;
    for k=1:M
        for v=1:N_TrialsLs2(k)
            Snoise = ST + NoiseIntensity*randn(N,M);
            Xnoise_UN(:,c) = Snoise(:,k);
            Onoise_UN(c) = k+Bias;
            c = c+1;
        end
    end
    
    XiNoisy_UN = (Onoise_UN*pinv(Xnoise_UN))';
    projs2(n,:) = ST'*XiNoisy_UN;
end
%%
color1 = [0 0.9 0.4];
color2 = 0.5*color1;
figure
plot(N_TrialsLs,'o-', LineWidth=2.5, MarkerSize=8, MarkerFaceColor='white', Color=color1)
hold on
plot(N_TrialsLs2,'o-', LineWidth=2.5, MarkerSize=8, MarkerFaceColor='white', Color=color2)
xticklabels(voc)
xlim([1 M])
ylabel("Training Trials")
pbaspect([1 1.5 1])
ax = gca;
set(ax, 'TickDir', 'out')
pbaspect([1 1.5 1])
if SAVE_FIGURES
    hgexport(gcf,'XiDifNoise_Trials1');
end
%%
figure
%plot(projs2', '-', LineWidth=2.5, MarkerSize=8, MarkerFaceColor='white', Color=[0.25*color2, 0.03])
hold on
%plot(projs', '-', LineWidth=2.5, MarkerSize=8, MarkerFaceColor='white', Color=[0.25*color1, 0.03])
%plot(mean(projs2), 'o-', LineWidth=2.5, MarkerSize=6, MarkerFaceColor='white', Color=color2)
%plot(mean(projs), 'o-', LineWidth=2.5, MarkerSize=6, MarkerFaceColor='white', Color=color1)
errorbar(1:M,mean(projs2), std(projs2)/sqrt(NNET), 'o-', LineWidth=2.5, MarkerSize=6, MarkerFaceColor='white', Color=color2)
errorbar(1:M,mean(projs), std(projs)/sqrt(NNET), 'o-', LineWidth=2.5, MarkerSize=6, MarkerFaceColor='white', Color=color1)
xticklabels(voc)
xlim([1 M])
%ylim([0 M-2])
grid off
ylabel('Projection')
title("<S|XiNoisy-ti>")
ax = gca;
set(ax, 'TickDir', 'out')
pbaspect([1 1.5 1])
if SAVE_FIGURES
    hgexport(gcf,'XiDifNoise_projs21');
end