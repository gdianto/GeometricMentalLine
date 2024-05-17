function [proj,Acc,AccS] = GML_ProbRewards(N,M,sigma,N_TrialsL,N_TrialsT,noiseL,noiseT)
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

% Orthogonalizing
[ST,~] = mgson(S);

% Adding Noise on couples

Xnoise = nan(N,N_TrialsL*(2*(M-1)));
Onoise = ones(1,size(Xnoise,2));
c = 1;
for k=1:(M-1)
    for v=1:N_TrialsL
        Snoise = ST + noiseL*randn(N,M);
        vec = Snoise(:,k+1)-Snoise(:,k);
        Xnoise(:,c) = vec;
        %Onoise(c) = log(1+1/k)+sigma*randn();
        Onoise(c) = sign(log(1+1/k)+sigma*randn());
        c = c+1;
    end
end
Xnoise(:,N_TrialsL*(M-1)+1:end) = -Xnoise(:,1:N_TrialsL*(M-1));

Onoise(N_TrialsL*(M-1)+1:end) = -Onoise(1:N_TrialsL*(M-1));

XiNoisy_ti = (Onoise*pinv(Xnoise))';

proj = ST'*XiNoisy_ti;
% Test on all couples
XnoiseSD = nan(N,M*(M-1)*N_TrialsT);
OnoiseSD = nan(1,M*(M-1)*N_TrialsT);
ImageNoiseR = nan(1,M*(M-1)*N_TrialsT);
ImageNoiseL = nan(1,M*(M-1)*N_TrialsT);
counter = 1;
for d=1:M
    for s=1:M
        sd = (s-d); %sdist
        if(sd~=0)
            for v=1:N_TrialsT
                Snoise = ST + noiseT*randn(N,M);
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

% Accuracy check
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

end