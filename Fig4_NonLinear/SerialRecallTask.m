function [projs,AccItem,AccSequence] = SerialRecallTask(M,N,stdJ,N_TrialsL,N_TrialsT,noise0)

S = randn(N,M);

%normalize vectors
for k=1:M
    S(:,k) = S(:,k)/norm(S(:,k));
end

% Orthogonalizing
[ST,~] = mgson(S);

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
projs = ST'*XiNoisyDelayed_single;

% Test memory and noise on single item
XSnoiseDelayed_T = nan(N,M*N_TrialsT);
for v=1:N_TrialsT
    Snoise = zeros(size(ST));
    for k=1:M
        Snoise(:,k) = (ST(:,k) + (noise0+(stdJ*(k-1)))*randn(N,1))/(1+(noise0*stdJ*k));
    end
    XSnoiseDelayed_T(:,(v-1)*M+1:M*v) = Snoise;
end

OSnoise_T = repmat([1:M],1,N_TrialsT);
% Single item performance
Resp = XSnoiseDelayed_T'*XiNoisyDelayed_single;
AccItem = zeros(1,M);
ndx = find(OSnoise_T == 1);
AccItem(1) = sum(Resp(ndx) < (1+0.5)) / numel(ndx);
ndx = find(OSnoise_T == M);
AccItem(M) = sum(Resp(ndx) > (M-0.5)) / numel(ndx);
for k=2:M-1
    ndx = find(OSnoise_T == k);
    AccItem(k) = sum((Resp(ndx) < (k+0.5)) & (Resp(ndx) > (k-0.5))) / numel(ndx);
end

% Global sequence performance
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
AccSequence = mean(correct_trials);

end