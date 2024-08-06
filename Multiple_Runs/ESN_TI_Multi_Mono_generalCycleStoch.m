clear all
close all
%%
M = 50; % number of ESNs to be tested

SEEDS = randi(100000, [1, M]);

Options.Tau = 0.1; %seconds
Options.EqualizeJinter = 0;
Options.NoiseForLearning = 0.05;

InputOptions.Wait = 30*Options.Tau; % wait time from one couple to the next one
InputOptions.Emission = 15*Options.Tau; % input signal width
InputOptions.EmissionR = 15*Options.Tau; % response signal width (+1=right,-1=left)
InputOptions.Cycles = 10000; % tMax ~ (wait+emission)*cycles
InputOptions.N_AP = 25;
InputOptions.N_SP = 30/7;

Options.Discrete = true;
Options.dt = 0.1 * Options.Tau;

% list of parameters to be tested
Ns = [100];
nInps = [7];
N_ALs = [20];
gs = [0.05];
ratio_g_gins = [5];
ratio_gin_gbacks = [999];
DELAYS = [0.1]*Options.Tau;
DNIs = 1.2:0.05:1.3;
SNIs = 0.0:0.05:1.3;
THs = [0.5];



fileName = ['TICycle(',...
    'N' num2str(Ns(1)) 'to' num2str(Ns(end)),...
    '_NINPS' num2str(nInps(1)) 'to' num2str(nInps(end)),...
    '_N_AL' num2str(N_ALs(1)) 'to' num2str(N_ALs(end)),...
    '_g' num2str(gs(1)) 'to' num2str(gs(end)),...
    '_rGGin' num2str(ratio_g_gins(1)) 'to' num2str(ratio_g_gins(end)),...
    '_rGinGb' num2str(ratio_gin_gbacks(1)) 'to' num2str(ratio_gin_gbacks(end)),...
    '_Delay' num2str(DELAYS(1)) 'to' num2str(DELAYS(end)),...
    '_DNI' num2str(DNIs(1)) 'to' num2str(DNIs(end)),...
    '_SNI' num2str(SNIs(1)) 'to' num2str(SNIs(end)),...
    '_TH' num2str(THs(1)) 'to' num2str(THs(end)),...
    ').csv'];

%--------------------

VarNames = ["nInp", "g", "rGGin", "rGinGb", "Delay", "N_AL", "DNI", "SNI","TH", "N"];
VarCells = cell(1,numel(VarNames));
VarCells{1} = nInps;
VarCells{2} = gs;
VarCells{3} = ratio_g_gins;
VarCells{4} = ratio_gin_gbacks;
VarCells{5} = DELAYS;
VarCells{6} = N_ALs;
VarCells{7} = DNIs;
VarCells{8} = SNIs;
VarCells{9} = THs;
VarCells{10} = Ns;


seqFileName = ['Seq_', fileName];
fileSeq = fopen(seqFileName,'w');
%fprintf(fileSeq, [columns, '\n']);
for k=1:numel(VarCells)
    obj = VarCells{k};
    fprintf(fileSeq, VarNames(k) + " ");
    for v=1:(numel(obj)-1)
        fprintf(fileSeq, '%g ', obj(v));
    end
    fprintf(fileSeq, '%g\n', obj(end));
end
fclose(fileSeq);


%-------------------------------
columns = "";
for name = VarNames
    columns = columns + name + " ";
end

for k=1:nInps(end)
    columns = columns +  "MeanANC" +  num2str(k) + " ";
    columns = columns +  "StdANC" + num2str(k) + " ";
end

columns = columns + "MeanA StdA SuccA ";

for k=1:(nInps(end) - 1)
    columns = columns + "MeanASD" + num2str(k) + " ";
    columns = columns + "StdASD" + num2str(k) + " ";
    columns = columns + "SuccASD" + num2str(k) + " ";
end


for k=1:(nInps(end) - 1)
    columns = columns +  "MaxOutSD" + num2str(k) + " ";
    columns = columns + "StdMO" + num2str(k) + " ";
end

for k=1:(nInps(end) - 1)
    columns = columns + "ReactTimesSD" + num2str(k) + " ";
    columns = columns + "StdRT" +  num2str(k) + " ";
end

for k=1:(nInps(end) - 1)
    columns = columns + "NeuroPercSD" + num2str(k) + " ";
    columns = columns + "StdNP" +  num2str(k) + " ";
end

for k=1:nInps(end)
    columns = columns + "Wrs" + num2str(k) + " ";
    columns = columns + "StdW" +  num2str(k) + " ";
end

for k=1:nInps(end)
    columns = columns + "WrsC" + num2str(k) + " ";
    columns = columns + "StdWC" +  num2str(k) + " ";
end

columns = columns + "PercError StdPE EVPC1 StdEV1 MBA StdMBA MBRt StdMBRt MBANC StdMBANC MBCorr StdMBCorr MBAll StdMBAll PercAnsw StdPA";
columns = columns + "\n";

%-----------------------------
%%
fileID = fopen(fileName,'w');
fprintf(fileID, columns);

counter = 0;
cycles = numel(Ns)*numel(N_ALs)*numel(gs)*numel(ratio_g_gins)*numel(ratio_gin_gbacks)*numel(DELAYS)*numel(SNIs)*numel(DNIs)*numel(nInps)*numel(THs);

ACC_TH = 0.9;
tic

for N=Ns
    for nInp = nInps
        for N_AL = N_ALs
            for TH = THs
                for DNI = DNIs
                    for SNI = SNIs
                        for g = gs
                            for r_ggin = ratio_g_gins
                                for r_gingb = ratio_gin_gbacks
                                    for DELAY=DELAYS
                                        Options.NetSize = N;
                                        Options.StdJinter = g;
                                        Options.StdJinput = g/r_ggin;
                                        if(r_gingb > 99)
                                            Options.StdJback = 0;
                                        else
                                            Options.StdJback = Options.StdJinput/r_gingb;
                                        end
                                        Options.DirNoise.Intensity = DNI * Options.StdJinput;
                                        Options.NoiseSD = SNI;
                                        Options.th = TH;
    
                                        InputOptions.NInputs = nInp;
                                        InputOptions.Delay = DELAY;
                                        InputOptions.N_AL = N_AL / nInp;
    
                                        fprintf(fileID, '%d %g %g %g %g %d %g %g %g %d ', nInp, g, Options.StdJinput,Options.StdJback, DELAY, N_AL, DNI, SNI, TH, N);
            
                                        %------------------------------
                                        counter = counter + 1;
                                        Anchors = zeros(M,nInps(end));
                                        AccSDs = zeros(M,nInps(end)-1);                           
                                        Accs = zeros(M,1);
                                        MaxOutSDs = zeros(M,nInps(end)-1);
                                        ReactTimesSDs = zeros(M,nInps(end)-1);
                                        NeuroPercSDs = zeros(M,nInps(end)-1);
                                        EVPC1 = zeros(M,1);
                                        ErrPjs = zeros(M,1);
                                        AccMB = zeros(M,1);
                                        RtMB = zeros(M,1);
                                        AnchMB = zeros(M,1);
                                        CorrMB = zeros(M,1);
                                        AllMB = zeros(M,1);
                                        Wrs = zeros(M,nInps(end));
                                        WrsC = zeros(M,nInps(end));
                                        PercAnswM = zeros(M,1);
            
                                        parfor v=1:numel(SEEDS)
                                            rng(SEEDS(v));
                                            fprintf('iteration = %d/%d-%d/%d ...\n', counter, cycles, v, M);
                                            [AnchorAccTemp, SymDist, MonkEffect, PCheck, Distr] = ESN_TI_Multi_MonoF_Stoch(Options, InputOptions);
                                            if(numel(AnchorAccTemp) < nInps(end)) % because of cycles on #inputs
                                                AnchorAcc = nan(1,nInps(end));
                                                SD_Acc = nan(1,nInps(end)-1);
                                                SD_MaxOut = nan(1,nInps(end)-1);
                                                SD_ReactTimes = nan(1,nInps(end)-1);
                                                SD_NeuroPerc = nan(1,nInps(end)-1);
                                                WrsTest = nan(1,nInps(end));
                                                WrsCTest = nan(1,nInps(end));
            
                                                AnchorAcc(1:numel(AnchorAccTemp)) = AnchorAccTemp;
                                                SD_Acc(1:numel(AnchorAccTemp)-1) = SymDist.Acc;
                                                SD_MaxOut(1:numel(AnchorAccTemp)-1) = SymDist.MaxOut;
                                                SD_ReactTimes(1:numel(AnchorAccTemp)-1) = SymDist.ReactTimes;
                                                SD_NeuroPerc(1:numel(AnchorAccTemp)-1) = SymDist.NeuroPercs;
                                                WrsTest(1:numel(AnchorAccTemp)) = Distr.wrs_im;
                                                WrsCTest(1:numel(AnchorAccTemp)) = Distr.wrsC_im;
            
                                            else
                                                AnchorAcc = AnchorAccTemp;
                                                SD_Acc = SymDist.Acc;
                                                SD_MaxOut = SymDist.MaxOut;
                                                SD_ReactTimes = SymDist.ReactTimes;
                                                SD_NeuroPerc = SymDist.NeuroPercs;
                                                WrsTest = Distr.wrs_im;
                                                WrsCTest = Distr.wrsC_im;
                                            end
            
                                            Anchors(v,:) =  AnchorAcc;
                                            AccSDs(v,:) = SD_Acc;
                                            Accs(v) = mean(SymDist.Acc);
                                            MaxOutSDs(v,:) = SD_MaxOut;
                                            ReactTimesSDs(v,:) = SD_ReactTimes;
                                            NeuroPercSDs(v,:) = SD_NeuroPerc;
            
                                            EVPC1(v) = PCheck.EV_PC1;
                                            ErrPjs(v) = PCheck.Err_Pj;
            
                                            AccMB(v) = MonkEffect.AccBool;
                                            RtMB(v) = MonkEffect.RtBool;
                                            AnchMB(v) = MonkEffect.AncBool;
                                            CorrMB(v) = MonkEffect.CorrBool;
                                            AllMB(v) = MonkEffect.All;
            
                                            Wrs(v,:) = WrsTest;
                                            WrsC(v,:) = WrsCTest;
                                            
                                            PercAnswM(v) = sum(SymDist.OTHCounter) / sum(SymDist.Counter) ;
                                        end
            
                                        for k=1:nInps(end)
                                            MeanAnchors = mean(Anchors(:,k));
                                            StdANC = std(Anchors(:,k));
                                            fprintf(fileID, '%g %g ', MeanAnchors, StdANC);
                                        end
            
                                        MeanA = mean(Accs);
                                        StdA = std(Accs);
                                        SuccA = sum(uint16(Accs>ACC_TH))/M;
                                        fprintf(fileID, '%g %g %g ', MeanA, StdA, SuccA);
            
                                        MeanASD = mean(AccSDs);
                                        StdASD = std(AccSDs);
                                        SuccASD = nan(1,InputOptions.NInputs-1);
                                        for k=1:(nInps(end)-1)
                                            SuccASD(k) = sum(uint16(AccSDs(:,k)>ACC_TH))/M;
                                            fprintf(fileID, '%g %g %g ', MeanASD(k), StdASD(k), SuccASD(k));
                                        end
            
                                        for k=1:(nInps(end)-1)
                                            fprintf(fileID, '%g %g ', mean(MaxOutSDs(:,k)), std(MaxOutSDs(:,k)));
                                        end
            
                                        for k=1:(nInps(end)-1)
                                            fprintf(fileID, '%g %g ', mean(ReactTimesSDs(:,k)), std(ReactTimesSDs(:,k)));
                                        end
            
                                        for k=1:(nInps(end)-1)
                                            fprintf(fileID, '%g %g ', mean(NeuroPercSDs(:,k)), std(NeuroPercSDs(:,k)));
                                        end
            
                                        for k=1:nInps(end)
                                            fprintf(fileID, '%g %g ', mean(Wrs(:,k)), std(Wrs(:,k)));
                                        end
            
                                        for k=1:nInps(end)
                                            fprintf(fileID, '%g %g ', mean(WrsC(:,k)), std(WrsC(:,k)));
                                        end
            
                                        PercepError = mean(ErrPjs);
                                        StdPE = std(ErrPjs);
            
                                        EV1 = mean(EVPC1);
                                        StdEV1 = std(EVPC1);
            
                                        MBAcc = mean(AccMB);
                                        StdMBAcc = std(AccMB);
                                        MBRt = mean(RtMB);
                                        StdMBRt = std(RtMB);
                                        MBAnch = mean(AnchMB);
                                        StdMBAnch = std(AnchMB);
                                        MBCorr = mean(CorrMB);
                                        StdMBCorr = std(CorrMB);
                                        MBAll = mean(AllMB);
                                        StdMBAll = std(AllMB);
                                        PercAnsw = mean(PercAnswM);
                                        StdPA = std(PercAnswM);
            
                                        fprintf(fileID, '%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n',...
                                                PercepError, StdPE,...
                                                EV1, StdEV1,...
                                                MBAcc, StdMBAcc, MBRt, StdMBRt, MBAnch, StdMBAnch, MBCorr, StdMBCorr, MBAll, StdMBAll, PercAnsw, StdPA);
            
                                        %------------------------------
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
toc

fclose(fileID);
