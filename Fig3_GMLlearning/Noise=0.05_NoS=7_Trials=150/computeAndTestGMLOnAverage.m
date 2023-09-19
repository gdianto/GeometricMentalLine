%%  Code simulating a Perceptron-like shallow network learning to perform 
%   a transitive inference task. Figures resulting from this simulation and
%   the related analysis provided the material to compose Fig. 3 of 
%   (Di Antonio et al., bioRxiv 2023, DOI: 10.1101/2023.08.03.551859).
%
%   Copyright Maurizio Mattia - Mar. 20, 2023
%

%% Set params and Symbols.
NoS = 7;   % Number of Symbols.
NoU = 1000; % Number of units in the network.
NoT = 150; % Number of trials.
NoiseDir = 0.05; % Directional noise.
NoR = 20; % Numerb of repetitions.

Rank = (NoS:-1:1)';
Symbols = randn([NoU NoS]);
Labels = char('A'+(0:NoS-1));
for k = 1:NoS
   Symbols(:,k) = Symbols(:,k)/norm(Symbols(:,k));
end
Biases = ones(1,NoS);

%%
allSD = zeros(NoS,NoS,NoR);
allAccuracy = zeros(NoS,NoS,NoR);
allGMLVsGMLOut = zeros(NoR,1);
allGMLVsGMLIn = zeros(NoR,1);

for nr = 1:NoR
   fprintf('Repetition n. %d of %d\n',nr,NoR)

   %% Make the consolidation-like learning dataset.
   Responses = 2*(randi(2,[NoT 1]) - 1.5);% Left = -1, Right = +1
   Pairs = randi(NoS-1,[NoT 1]);
   Pairs = [Pairs Pairs+1];
   RankDiff = Rank(Pairs(:,2)) - Rank(Pairs(:,1));

   for k = find(Responses ~= RankDiff)
      Pairs(k,:) = fliplr(Pairs(k,:));
   end

   X = sqrt(2)*NoiseDir*randn(NoU,NoT);
   X = [X; ones(1,NoT)];
   for k = 1:NoT
      X(1:NoU,k) = -Symbols(:,Pairs(k,1))+Symbols(:,Pairs(k,2))+X(1:NoU,k);
   end

   GML = pinv(X)'*Responses; % Perceptron weights, i.e., the geometric mental line (GML).

   %% Compute the projections of pair of Symbols onto the GML, i.e., the
   %  signed symbolic distance (SD).
   x = ones(NoU+1,1);

   for l = 1:NoS
      for r = 1:NoS
         x(1:NoU,:) = Symbols(:,r) - Symbols(:,l);
         allSD(l,r,nr) = GML'*x;
      end
   end

   %% Compute the cosine of GML and GMLout
   Prjs = GML'*[Symbols; Biases];
   GMLin = [Symbols; Biases]*Prjs';
   GMLout = GML-GMLin;

   allGMLVsGMLIn(nr) = GMLin'*GML/(norm(GMLin)*norm(GML));
   allGMLVsGMLOut(nr) = GMLout'*GML/(norm(GMLout)*norm(GML));

   %% Estimate and plot accuracy per symbol pair.
   NoAT = 100000; % Num. of trials to test estimate accuracy

   BGNoise = (GML(1:NoU)'*(sqrt(2)*NoiseDir*randn(NoU,NoAT)));
   x = ones(NoU+1,1);
   for l = 1:NoS
      for r = 1:NoS
         if l ~= r
            x(1:NoU,:) = Symbols(:,r) - Symbols(:,l);
            y = GML'*x+BGNoise;
            y = (Rank(r)-Rank(l))*y;
            allAccuracy(l,r,nr) = numel(find(y>0))/NoAT;
         end
      end
   end

end % for nr = ...

SD = mean(allSD,3);
Accuracy = mean(allAccuracy,3);

%% Plot the SD matrix and its relationship with the projections on the GML.
figure
subplot(1,2,1)
imagesc(SD,[-1 1]*NoS)
colormap(RedBlueCM)
% title(sprintf('Estimated SD, N = %d, Trials = %d, Noise = %g',NoU,NoT,NoiseDir))
title('Estimated SD')
xlabel('Right symbol')
ylabel('Left symbol')
colorbar

hold on

for l = 1:NoS
   for r = 1:NoS
      x(1:NoU,:) = Symbols(:,r) - Symbols(:,l);
      text(r,l,num2str(SD(l,r),'%.2f'),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',9);
   end
end

set(gca,'DataAspectRatio',[1 1 1],'TickDir','out')
xticks(1:NoS)
xticklabels(Labels')
yticks(1:NoS)
yticklabels(Labels')

% Plot the GML projections versus SD.
%figure
subplot(1,2,2)
plot([-1 1]*(NoS-1),[-1 1]*(NoS-1),'k-')
hold on
x = zeros(NoS*(NoS-1),1);
y = x;
k = 0;
for r = 1:NoS
   for l = 1:NoS
     if r ~= l
        k = k + 1;
        x(k) = Rank(r)-Rank(l);
        y(k) = SD(l,r);
     end
   end
end
plot(x,y,'ko','MarkerFaceColor','w')
pl = polyfit(x,y,1);
plot([-1 1]*(NoS-1),polyval(pl,[-1 1]*(NoS-1)),'r--')

xlabel('Symbolic Distance (SD)')
ylabel('GML projection')
title(sprintf('N = %d, Trials = %d, Noise = %g, Slope = %.2g',NoU,NoT,NoiseDir,pl(1)))

set(gca,'DataAspectRatio',[1 1 1],'TickDir','out')

FigSize = [10 3.5];
set(gcf, 'PaperUnits', 'inch', 'PaperSize', FigSize, 'PaperPosition', [0 0 FigSize]);
print('-dpdf', 'averageSDandGMLprojections');

%% Plot the Accuracy matrix.
figure

subplot(2,2,3)
imagesc(Accuracy,[0 1])
colormap(bone)
% title(sprintf('Accuracy, N = %d, Trials = %d, Noise = %g',NoU,NoT,NoiseDir))
title('Accuracy')
xlabel('Right symbol')
ylabel('Left symbol')
colorbar

hold on

for l = 1:NoS
   for r = 1:NoS
      x(1:NoU,:) = Symbols(:,r) - Symbols(:,l);
      text(r,l,num2str(Accuracy(l,r),'%.2f'),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',9);
   end
end

set(gca,'DataAspectRatio',[1 1 1],'TickDir','out')
xticks(1:NoS)
xticklabels(Labels')
yticks(1:NoS)
yticklabels(Labels')

% PlotAccuracy versus Symbol pairs.
Z = (Accuracy+Accuracy')/2;

subplot(2,2,1:2)
hold on
Offset = 0;
XTicks = [];
XTickLabels = [];
for j = 1:NoS-1 % Loop on SD.
   x = (1:NoS-j) + Offset;
   y = 0*x;
   for k = 1:NoS-j
      y(k) = Z(k,k+j);
      XTickLabels = [XTickLabels; Labels([k k+j])];
   end
   plot(x,y,'k-d','LineWidth',0.75)
   Offset = x(end)+1;
   XTicks = [XTicks x];
end
xlim([0 Offset+1])
ylim([0 1.1])
xticks(XTicks)
xticklabels(XTickLabels)
title(sprintf('N = %d, Trials = %d, Noise = %g',NoU,NoT,NoiseDir))

ylabel('Accuracy')
set(gca,'TickDir','out')
grid on

% Mean accuracy per symbol.
subplot(2,2,4)
plot(1:NoS,sum(Accuracy)/(NoS-1),'k-d')
ylim([0 1.1])
xticks(1:NoS)
xlim([0 NoS+1])
xticklabels(Labels)
xticklabels(Labels')
grid on

FigSize = [8 6];
set(gcf, 'PaperUnits', 'inch', 'PaperSize', FigSize, 'PaperPosition', [0 0 FigSize]);
print('-dpdf', 'averageAccuracyAndSD');

%% Plot Accuracy versus SD.
figure
%subplot(1,2,2)
hold on
x = zeros(NoS*(NoS-1),1);
y = x;
k = 0;
for r = 1:NoS
   for l = 1:NoS
     if r ~= l
        k = k + 1;
        x(k) = abs(Rank(r)-Rank(l));
        y(k) = Accuracy(l,r);
     end
   end
end
plot(x,y,'ko','MarkerFaceColor','w')

xm = 1:max(x);
ym = zeros(size(xm));
for nSD = 1:max(x)
   ym(nSD) = mean(y(x==nSD));
end
SDrange = [1 NoS-1];
% pl = polyfit(x,y,3);
% plot(linspace(SDrange(1),SDrange(2),(NoS-1)*4),polyval(pl,linspace(SDrange(1),SDrange(2),(NoS-1)*4)),'r--')
plot(xm,ym,'ro--','LineWidth',0.75,'MarkerFaceColor','w')

xlim(SDrange+[-1 1]*0.1)
ylim([0 1.05])
grid on

xlabel('|SD|')
ylabel('Accuracy')
title(sprintf('N = %d, Trials = %d, Noise = %g',NoU,NoT,NoiseDir))
set(gca,'TickDir','out')

FigSize = [4 4];
set(gcf, 'PaperUnits', 'inch', 'PaperSize', FigSize, 'PaperPosition', [0 0 FigSize]);
print('-dpdf', 'averageAccuracyVsAbsSD');

%% Plot boxplot with cosines of GML and GMLout.
figure
% boxplot(allGMLVsGMLOut,'Notch','on')
boxplot(allGMLVsGMLIn,'Notch','on')
% ylabel('\langle\zeta|\zeta_{out}\rangle')
ylabel('\langle\zeta|\zeta_{in}\rangle')
title(sprintf('N = %d, Trials = %d, Noise = %g',NoU,NoT,NoiseDir))
ylim([0 1])
grid on
set(gca,'TickDir','out')
xticklabels({''})
xlabel(sprintf('\\sigma_D = %.3g',NoiseDir))

FigSize = [4 4];
set(gcf, 'PaperUnits', 'inch', 'PaperSize', FigSize, 'PaperPosition', [0 0 FigSize]);
% print('-dpdf', 'averageGMLxGMLout');
print('-dpdf', 'averageGMLxGMLin');
