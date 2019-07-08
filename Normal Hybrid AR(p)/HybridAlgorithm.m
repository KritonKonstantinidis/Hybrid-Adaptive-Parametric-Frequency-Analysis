function [] = HybridAlgorithm(y,varargin)

% This function takes as input a raw time-series y and produces its 
% time-varying spectrogram using autoregressive models, online EM algorithm 
% and hybrid Kalman filtering.

%% Default values of parameters
opt.exp_smoothing=0; % exponential smoothing of AR coefficients
opt.fs = 250; % sampling rate 
opt.psdresolution = 200; % frequency resolution 
opt.maxfreq=50; % maximum of the frequency range
opt.minfreq=0.00001; % minimum of the frequency range
opt.maxIter=50; % maximum number of EM algorithm
opt.minibatch_size=1250; %  5 seconds of EEG data
opt.tol_EM=1e-3; % tolerance for EM convergence
opt.Rstart=0.15; % initial observation noise variance for EM iterations
opt.Qstartcoef=1; % starting coefficient of diagonal covariance matrix for EM iterations
opt.pmin=2; % minimum AR order to be considered
opt.pmax=20; % maximum AR order to be considered
opt.criterion= 1; % information criterion (0 for BIC or 1 for AIC)
opt.remove=1; % remove outliers or not 
opt.online=0; % update Q or not
opt.buffer_size=30000; % buffer size in case of online EM
%% Parse input arguments 
opt = parsevarargin(varargin, opt); 

exp_smoothing=opt.exp_smoothing;
fs=opt.fs;
psdresolution=opt.psdresolution;
maxfreq=opt.maxfreq;
minfreq=opt.minfreq;
maxIter=opt.maxIter;
minibatch_size=opt.minibatch_size;
Rstart=opt.Rstart;
Qstartcoef=opt.Qstartcoef;
pmin=opt.pmin;
pmax=opt.pmax;
criterion=opt.criterion;
tol_EM=opt.tol_EM;
remove=opt.remove;
online=opt.online;
buffer_size=opt.buffer_size;
%% Check the inputs
if exp_smoothing ~=0 && exp_smoothing ~=1
    error('Please input exponetnial smoothing (1) or not (0)');
end
if criterion ~=0 && criterion ~=1
    error('Please input BIC (0) or AIC (1)');
end
if online ~=0 && online ~=1
    error('Please input adaptive estimation of Q (1) or not (0)');
end
if remove ~=0 && remove ~=1
    error('Please input outlier rejection (1) or not (0)');
end
%% Remove Outliers
if remove==1
N=size(y,1);
threshold=mean(y)+5*std(y);
for i=1:N
    if abs(y(i))>threshold
        y(i)=0;
    end
end
y(y==0) = [];
end
%% Sample Number and Data duration
N=size(y,1); % number of samples
duration=N/fs; % duration of the timeseries in seconds
nsec=duration/5; % show x axis ticks/plot every duration/5 seconds
%% Data Normalization
maxY = max(abs(y));
y = y./maxY;
%% Optimization of noise and autoregressive model order
[p,Q,R,P0] = ModelSelection(pmin,pmax,criterion,y,minibatch_size,Rstart,Qstartcoef,fs,exp_smoothing,maxIter,tol_EM);

fprintf(['Selected order of the autoregressive model ',num2str(p)]);
fprintf('\n') ;
fprintf('EM Process Noise Covariance : ') ; Q
fprintf('\n') ;
fprintf(['EM Observation Noise Variance:',' ', num2str(R)]);
%% Hybrid Kalman Filtering
if online==0
    [afilt,apreds,Pfilt,Ppreds,S]=HybridKalmanFilter(y,p,P0,R,Q,fs,exp_smoothing);
    a=HybridKalmanSmoother(y,p,afilt,apreds,Pfilt,Ppreds,S);
end

if online==1
    Qstart=Qstartcoef*eye(p);
    ainit=aryule(y,p);
    aprev=-ainit(2:end)';
    afilt=cell(ceil(N/buffer_size),1);
    apreds=cell(ceil(N/buffer_size),1);
    Pfilt=cell(ceil(N/buffer_size),1);
    Ppreds=cell(ceil(N/buffer_size),1);
    Q=cell(ceil(N/buffer_size),1); 
    S=cell(ceil(N/buffer_size),1); 

    %%% 1st pass of the algorithm %%%%
    [~,Q{1},~]=HybridEMwithQ(y(1:minibatch_size),p,P0,maxIter,R,Qstart,fs,exp_smoothing,tol_EM);
    [afilt{1},apreds{1},Pfilt{1},Ppreds{1},S{1}]=HybridKalmanFilterOnline(y(1:buffer_size),p,P0,R,Q{1},fs,exp_smoothing,aprev);
    aprev=afilt{1}(:,end);
    j=2;

    for i=buffer_size+1:buffer_size:N
        if i+minibatch_size>N
            minibatch_size=N-i;
        end
        if i+buffer_size>N
            buffer_size=N-i+1;
        end

        %%% EM in minibatch
        [~,Q{j},~]=HybridEMwithQ(y(i:i+minibatch_size),p,P0,maxIter,R,Qstart,fs,exp_smoothing,tol_EM);
        %%% Run forward Kalman 
        [afilt{j},apreds{j},Pfilt{j},Ppreds{j},S{j}]=HybridKalmanFilterOnline(y(i:i+buffer_size-1),p,P0,R,Q{j},fs,exp_smoothing,aprev);
        aprev=afilt{j}(:,end);
        j=j+1;
    end

    afilt=horzcat(afilt{:});
    apreds=horzcat(apreds{:});
    Pfilt=cat(3,Pfilt{:});
    Ppreds=cat(3,Ppreds{:});
    S=cat(3,S{:});
    a=HybridKalmanSmoother(y,p,afilt,apreds,Pfilt,Ppreds,S);
end
%% Spectrogram
specgram = HybridSpectrogram(a,minfreq, maxfreq, R,'psdres',psdresolution); 
figure
imagesc(specgram.timestamp, specgram.w, specgram.specgram);
xlim([0 size(y,1)])
set(gca,'XTick',1:nsec*fs:N); 
set(gca,'XTickLabel', 0:nsec/60:N/fs/60);
xlabel('Time (Minutes)');
caxis([-15 10]);
colorMap=jet;
colormap(colorMap);
colorBar=colorbar;
ylabel(colorBar,'Power/Frequency (dB/Hz)');
axis xy;
ylabel('Freq (Hz)');
title('Hybrid Kalman Filter Spectrogram');
%% Plot coefficient evolution
a=a';
figure
for i=1:p
    subplot(ceil(p/2),2,i)
    plot(a(:,i))
    xlim([0 size(y,1)])
    set(gca,'XTick',1:nsec*fs:N); 
    set(gca,'XTickLabel', 0:nsec:N/fs);
    xlabel('Time (seconds)');
    ylabel('Amplitude');
    title(['a',num2str(i),' ','continuous estimation']);
end
end