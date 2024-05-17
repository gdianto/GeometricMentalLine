%%  Code simulating a the stochastic dynamics of the projection on the 
%   geometric mental line of the network state receiving as input the 
%   symbolic distance of the presente pair of items.
%   Figures resulting from this simulation/script and the related analysis 
%   provided the material to compose Fig. 4 of 
%   (Di Antonio et al., bioRxiv 2023, DOI: 10.1101/2023.08.03.551859).
%
%   Copyright Maurizio Mattia - Mar. 10, 2023
%

%% Params and colors.
% NoiseSize = 0.075*sqrt(2);
% NoiseSize = 0.5*sqrt(2);
NoiseSize = 0.75;
RandomSamples = 5;

Example = 1;
switch Example
   case 1
      SymbolRank = [4 3 1]; % A B D
      SDCM = SymbDistCM(max(SymbolRank));
end

%% PDF of OU process
RelPeakSize = 0.15;
NormPDF = @(x,m,s) RelPeakSize*exp(-(x-m).^2/(2*s^2))/(s*sqrt(2*pi));

Tau = 0.1;
Life = 0.35;

Mu = @(t,m,tau) m*(1-exp(-t/tau));
Sigma = @(t,s,tau) sqrt(s^2*(1-exp(-(2*t)/tau)));

NormPrctl = norminv([0.075 0.2:0.1:0.8 0.925])';
NoP = numel(NormPrctl);

figure
hold on

for SD = [1 -3]
   MInf = SD;
   SInf = NoiseSize;
   t = linspace(0,Life,100);
   y = zeros(numel(NormPrctl),numel(t));
   for k = 1:numel(t)
      y(:,k) = NormPrctl * Sigma(t(k),SInf,Tau) + Mu(t(k),MInf,Tau);
   end

   MaxK = floor(NoP/2);
   for k = 1:MaxK
      Color = SDCM(SD + max(SymbolRank) + 1,:);
      Color = [1 1 1]*(MaxK-k+3)/(MaxK+3) + Color*k/(MaxK+3);
      patch([t fliplr(t)],[y(k,:) fliplr(y(NoP-k+1,:))],Color,'EdgeColor','none')
   end
   plot(t,y(5,:),'-','LineWidth',0.75,'Color',SDCM(SD + max(SymbolRank) + 1,:))

   CurvAbsc = linspace(MInf-3.5*SInf,MInf+3.5*SInf,100);
   plot(Life - NormPDF(CurvAbsc,MInf,SInf),CurvAbsc,'-','LineWidth',0.75,'Color',SDCM(SD + max(SymbolRank) + 1,:))
end
set(gca,'TickDir','out')
xlabel('Time, t [s]')
ylabel('GML integrator, y(t)')
grid on 
xlim([0 Life])
ylim([-6 4])

FigSize = [4 6];
set(gcf,'PaperUnit','inch','PaperPosition',[0 0 FigSize],'PaperSize',FigSize);
print('-dpdf', '-painters', 'IntegratorOUVsTime');

%% OU Integrator sample realizations.
Life = 0.75;
dt = 1/1000;
Threshold = 0.8;

t = 0:dt:Life;

figure
hold on
for SD = 0 % [-3 1]
   MInf = SD;
   SInf = NoiseSize;
   for n = 1:RandomSamples
      x = 0*t;
      GWN = randn(size(t))*sqrt(2*dt/Tau);
      for k = 1:numel(t)-1
         x(k+1) = ((-x(k) + MInf)*dt/Tau + SInf*GWN(k)) + x(k);
      end
      plot(t,x,'-')
   end
end
plot([0 Life],[1 1]*Threshold,'k--')
plot([0 Life],-[1 1]*Threshold,'k--')
set(gca,'TickDir','out')
xlabel('Time, t [s]')
ylabel('GML integrator, y(t)')
grid on 
xlim([0 Life])
ylim([-6 4])

FigSize = [4 6];
set(gcf,'PaperUnit','inch','PaperPosition',[0 0 FigSize],'PaperSize',FigSize);
print('-dpdf', '-painters', sprintf('IntegratorOUVsTime_SD=%d',SD));
