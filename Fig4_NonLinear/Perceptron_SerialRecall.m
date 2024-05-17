%% 
clear all
close all
%% Create mean vectors (Symbols)
SAVE_FIGURES = false;
M = 6; % #symbols
N = 100; % dimension
S = randn(N,M);

alp = ('A':'Z')';
%voc = alp(1:M);
voc = alp(M:-1:1);

%normalize vectors
for k=1:M
    S(:,k) = S(:,k)/norm(S(:,k));
end

%% Orthogonalizing
[ST,~] = mgson(S);

%% Adding memory and noise on single item
N_TrialsL = 10000; % trials per item
gen_mol = 1;
noise0 = gen_mol*0.02;
stdJ = gen_mol*0.0075;
%stdJ = 0.02;
%stdJ = 0.12;

XSnoiseDelayed = nan(N,M*N_TrialsL);
for v=1:N_TrialsL
    Snoise = zeros(size(ST));
    for k=1:M
        Snoise(:,k) = (ST(:,k) + (noise0+(stdJ*(k-1)))*randn(N,1))/(1+(noise0*stdJ*k));
    end
    XSnoiseDelayed(:,(v-1)*M+1:M*v) = Snoise;
end

OSnoise = repmat([1:M],1,N_TrialsL);
XiNoisyDelayed_single = (OSnoise*pinv(XSnoiseDelayed))';
%%
figure
scatter(OSnoise, XiNoisyDelayed_single' * XSnoiseDelayed, 300,'.')
ylabel("proj. on delayed S")
xticklabels(voc)
pbaspect([1 1.2 1])
grid on

%% Test memory and noise on single item
N_TrialsT = 10000; % trials per item
%noise0 = 0.03;

XSnoiseDelayed_T = nan(N,M*N_TrialsT);
for v=1:N_TrialsT
    Snoise = zeros(size(ST));
    for k=1:M
        Snoise(:,k) = (ST(:,k) + (noise0+(stdJ*(k-1)))*randn(N,1))/(1+(noise0*stdJ*k));
    end
    XSnoiseDelayed_T(:,(v-1)*M+1:M*v) = Snoise;
end

OSnoise_T = repmat([1:M],1,N_TrialsT);

%%
all_proj = XiNoisyDelayed_single' * XSnoiseDelayed_T;
Y = zeros(N_TrialsT,M);
for k=1:M
    Y(:,k) = all_proj(OSnoise_T == k)';
end

%%
colors = [255 51 51;
    250 201 56;
    168 245 61;
    45 211 86;
    139 243 243;
    76 114 230;
    153 81 225]/255;

alph = 0.3;
figure
%scatter(OSnoise_T, all_proj, 150,'k.')
violin(Y,'facecolor',colors,'edgecolor','k',...
'bw',0.3,...
'mc','k',...
'linewidth', 5,...
'plotlegend', false,...
'medc','',14)

hold on
%rectangle('Position', [1, -5.5, M, 7], 'FaceColor',[colors(1,:) alph])
plot([1-0.5, 1+0.5],[1+0.5,1+0.5], '--', 'Color',colors(1,:), 'LineWidth',2)
scatter(1, 1, 2000,'w.')
scatter(1, 1, 1000,'.', 'MarkerEdgeColor',colors(1,:))

%rectangle('Position', [1, M-0.5, M, 10], 'FaceColor',[colors(M,:) alph])
%yline(M-0.5, '--', 'Color',colors(M,:), 'LineWidth',2)
plot([M-0.5, M+0.5],[M-0.5,M-0.5], '--', 'Color',colors(M,:), 'LineWidth',2)
scatter(M, M, 2000,'w.')
scatter(M, M, 1000,'.', 'MarkerEdgeColor',colors(M,:))

for k=2:M-1
    %rectangle('Position', [1, k-0.5, M, 1], 'FaceColor', [colors(k,:) alph])
    plot([k-0.5, k+0.5],[k+0.5,k+0.5], '--', 'Color',colors(k,:), 'LineWidth',2)
    plot([k-0.5, k+0.5],[k-0.5,k-0.5], '--', 'Color',colors(k,:), 'LineWidth',2)
    scatter(k, k, 2000,'w.')
    scatter(k, k, 1000,'.', 'MarkerEdgeColor',colors(k,:))
end
xlim([1,M+0.05])
ylim([-3,12])
ylabel("projection")
xlabel('Serial Position')
%xticklabels(voc)
ax = gca;
set(ax, 'TickDir', 'out', 'Box', 'off')
pbaspect([1.0 1.5 1])
if SAVE_FIGURES
    %hgexport(gcf,"XiMem_proj_"+num2str(stdJ) + ".eps");
    print('-painters','-dpdf',"XiMem_proj_"+num2str(stdJ) + ".pdf")
end
%% Single item performance
Resp = XSnoiseDelayed_T'*XiNoisyDelayed_single;
Acc = zeros(1,M);
ndx = find(OSnoise_T == 1);
Acc(1) = sum(Resp(ndx) < (1+0.5)) / numel(ndx);
ndx = find(OSnoise_T == M);
Acc(M) = sum(Resp(ndx) > (M-0.5)) / numel(ndx);
for k=2:M-1
    ndx = find(OSnoise_T == k);
    Acc(k) = sum((Resp(ndx) < (k+0.5)) & (Resp(ndx) > (k-0.5))) / numel(ndx);
end
%% Plot single item performance
figure
plot(Acc,'o-k', LineWidth=2.5, MarkerSize=8, MarkerFaceColor='white')
ylim([0 1])
xlim([1 M])
xticklabels(voc)
grid on
ylabel('Accuracy')
ax = gca;
set(ax, 'TickDir', 'out')
pbaspect([1 1.5 1])
if SAVE_FIGURES
    hgexport(gcf,"XiMem_"+num2str(stdJ));
end
%% Global sequence performance (whole sequence correct)
correct_trials = zeros(1,N_TrialsT);
for v=1:N_TrialsT
    in_trial_resp = Resp((v-1)*M+1:M*v);
    bool_correct_trial = (in_trial_resp(1) < (1+0.5)) & (Resp(M) > (M-0.5));
    if bool_correct_trial
        for k=2:(M-1)
            bool_correct_trial = bool_correct_trial & (in_trial_resp(k) < (k+0.5)) & (in_trial_resp(k) > (k-0.5));
        end
    end
    correct_trials(v) = bool_correct_trial;
end
sequence_accuracy = mean(correct_trials)

%% Different slope
N = 100;
N_TrialsL=100;
N_TrialsT=100;
%noise0 = 0.0075;
noise0 = 0.0075;

M = 6;
REP = 100;
%stdJs = 0.01:0.01:0.04;
stdJs = [0.0075, 0.015, 0.03, 0.06, 0.12];
%stdJs = 0.01*[1, 1, 1];
%noises0 = noise0 * [1, 5, 10];
mean_acc_item = zeros(M,numel(stdJs));
std_acc_item = zeros(M,numel(stdJs));
mean_projs = zeros(M,numel(stdJs));
std_projs = zeros(M,numel(stdJs));
for k=1:numel(stdJs)
    accs = zeros(M,REP);
    projs = zeros(M,REP);
    for r=1:REP
        [projs(:,r),accs(:,r),~] = SerialRecallTask(M,N,stdJs(k),N_TrialsL,N_TrialsT,noise0);
        %[projs(:,r),accs(:,r),~] = SerialRecallTask(M,N,stdJs(k),N_TrialsL,N_TrialsT,noises0(k));
    end
    if REP > 1
        mean_acc_item(:,k) = mean(accs');
        std_acc_item(:,k) = std(accs');
        mean_projs(:,k) = mean(projs');
        std_projs(:,k) = std(projs');
    else
        mean_acc_item(:,k) = accs;
        mean_projs(:,k) = projs';
    end
end
%% plot
%hue = flip(brewermap(numel(stdJs)+5, 'Reds'));
hue = flip(brewermap(numel(stdJs)+5, 'Blues'));
hue = flip(hue(1:numel(stdJs)+3,:));

SAVE_FIGURES=false;

figure
hold on
for k=1:numel(stdJs)
    plot(1:M, (1-mean_acc_item(:,k))*100, "LineWidth", 3, "Color", hue(k,:))
    errorbar(1:M, (1-mean_acc_item(:,k))*100, std_acc_item(:,k)/sqrt(REP),"LineStyle","none", "LineWidth", 3,...
        "Marker", ".", "MarkerSize", 25, "Color", hue(k,:))
    plot(1:M, (1-mean_acc_item(:,k))*100, "LineStyle","none","Marker", ".", "MarkerSize", 12, "Color", 'w')
end
ylim([0 100])
xlim([1 M])
%xticklabels(voc)
grid off
xlabel("serial position")
ylabel("errors(%)")
pbaspect([1 1.5 1])

if SAVE_FIGURES
    hgexport(gcf,"SerialRecall_accItem_diffMemories");
end
%%
figure
hold on
for k=1:numel(stdJs)
    %errorbar(1:M,mean_projs(:,k), std_projs(:,k)/sqrt(REP), 'o-', LineWidth=2.5, MarkerSize=8, MarkerFaceColor='white', color=hue(k,:))
    plot(1:M,mean_projs(:,k), 'o-', LineWidth=2.5, MarkerSize=8, MarkerFaceColor='white', color=hue(k,:))
end
%xticklabels(voc)
xlim([1 M])
grid off
ylabel('Projection <GML|S>')
xlabel('Serial Position')
ax = gca;
set(ax, 'TickDir', 'out')
pbaspect([1 1.5 1])

if SAVE_FIGURES
    hgexport(gcf,"SerialRecall_projs_diffMemories.eps");
end
%% Different size and slope
N = 100;
N_TrialsL=100;
N_TrialsT=100;

Ms = 3:8;
REP = 100;
mean_acc_sequences = zeros(numel(Ms), numel(stdJs));
std_acc_sequences = zeros(numel(Ms), numel(stdJs));
acc_sequences = zeros(numel(Ms), REP);

for n=1:numel(stdJs)
    for k=1:numel(Ms)
        parfor r=1:REP
            [~,~,acc_sequences(k,r)] = SerialRecallTask(Ms(k),N,stdJs(n),N_TrialsL,N_TrialsT,noise0);
        end
    end
    if REP > 1
        mean_acc_sequences(:,n) = mean(acc_sequences');
        std_acc_sequences(:,n) = std(acc_sequences');
    else
        mean_acc_sequences(:,n) = acc_sequences;
    end
end
%%
hue = flip(brewermap(numel(stdJs)+5, 'Blues'));
hue = flip(hue(1:numel(stdJs)+3,:));

figure
hold on
for k=1:numel(stdJs)
    plot(Ms,mean_acc_sequences(:,k), "LineWidth", 3, "Color", hue(k,:))
    errorbar(Ms,mean_acc_sequences(:,k), std_acc_sequences(:,k)/sqrt(REP),"LineStyle","none", "LineWidth", 3,...
        "Marker", ".", "MarkerSize", 25, "Color", hue(k,:))
    plot(Ms,mean_acc_sequences(:,k), "LineStyle","none","Marker", ".", "MarkerSize", 12, "Color", 'w')
end
ylim([0 1.0])
xlim([Ms(1) Ms(end)])
grid off
xlabel("list length")
ylabel("proportion perfect recall")
pbaspect([1 1.2 1])

if SAVE_FIGURES
    hgexport(gcf,"SerialRecall_accSeq_diffMemories.eps");
end