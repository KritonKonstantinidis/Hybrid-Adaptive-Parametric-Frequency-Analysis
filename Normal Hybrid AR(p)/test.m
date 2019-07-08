%% Load of real EEG data
clear all
close all
clc

tic
fs = 250; % Sampling frequency [Hz] (index/sec)

% 
% %%%%% Sevoflurane data
% % 
% y=load(fullfile('Datasets','SED149.mat.mat'));
% y=y.data;
% y=y(2,:);
% dataEEG=y(124*60*250:174*60*250)';
% dataEEG=dataEEG(1:100000);
%%%% Ketamine data
% 
% y=load(fullfile('Datasets','SED974_lite.mat'));
% y=y.data;
% % % % %%% beta gamma oscillations
%  y=y(1,:); % channel 1 
% dataEEG=y(306*60*250:316*60*250)';
% dataEEG=dataEEG(1:1000);
% %%%% bursts 
%  y=y(3,:);
%  dataEEG=y(3200*250:3600*250-1)'; 
% dataEEG=y(200*250:600*250-1)'; 
% % % 
% % %%%% Propofol data 
% 
data = load(fullfile('Datasets','SED10.mat'));
data=data.data;
channel=1; % Frontal Channel
dataEEG = data(channel,600*fs+1:1880*fs)'; % alpha oscillations
% dataEEG=dataEEG(1:30000);

% % % % % 
% [S,freqs]=spectrogram(dataEEG,[],[],[],250,'yaxis');
% 
N=size(dataEEG,1); % number of observations
duration=N/fs; % duration in seconds
nsec=duration/5; % show x axis ticks/plot instantaneous spectra every duration/5 seconds
indexstep=nsec*fs; 
% % tic
% HybridAlgorithm(dataEEG,'pmin',2,'pmax',16,'maxIter',50,'minibatch_size',22000,'online',0,'Rstart',0.069,'maxfreq',30);
% toc

%%%%%% REMOVE OUTLIERS %%%%%%%
% threshold=mean(dataEEG)+5*std(dataEEG);
% for i=1:N
%     if abs(dataEEG(i))>threshold
%         dataEEG(i)=0;
%     end
% end
% dataEEG(dataEEG==0) = [];
% % 
% N=size(dataEEG,1);
%Normalization 
maxEEG = max(abs(dataEEG));
dataEEG = dataEEG./maxEEG;
% tic
% %% Estimate Q via EM on an initial EEG segment
% clc 
% close all
% exp_smoothing=0;
% minibatch=dataEEG(1:1250);
% tol=1e-3;
% nIter=100;
p=14;
P0=eye(p);
% Rstart=0.07;
% Qstart=eye(p);
% tic
% [R,Q,L]=DiscreteEM(minibatch,p,P0,nIter,Rstart,Qstart,exp_smoothing,tol);
% [R,Q,L]=HybridEMwithQ(minibatch,p,P0,nIter,Rstart,Qstart,fs,exp_smoothing,tol);
% toc
 %% Hybrid Kalman Filter by manually choosing the model order
%  R=0.07;% Variance of the discrete gaussian white observation noise
%  Q=7e-3*eye(p);% Covariance matrix of white process noise 
% % P0=eye(p);% Initial state covariance matrix 
% exp_smoothing=0; % exponential smoothing of the estimates 
% % % % % [a,~] = dnf_arkal(dataEEG,fs, 'ARorder', 12,'ksmooth',0);
% % % % 
% tic
% [a_discrete,~,ypred_discrete]=KalmanFilterDiscrete(dataEEG,p,P0,R,Q,0);
% toc
% %  tic
% %  F=-0.05*eye(p);
% tic
% [a_hybrid,~,~,~,~,~,ypred_hybrid]=HybridKalmanFilter(dataEEG,p,P0,R,Q,fs,exp_smoothing);
% toc
% % [a,apreds,Pfilt,Ppreds,S,E,ypred]=HybridKalmanFilter(dataEEG,p,P0,R,Q,fs,exp_smoothing);
% a=HybridKalmanFilter(dataEEG,p,P0,R,Q,fs,exp_smoothing);
% 
% %  a=HybridKalmanSmoother(dataEEG,p,a,apreds,Pfilt,Ppreds,S);
% % % 
% %% Plot model against raw 
% figure
% 
% subplot(2,1,1)
% plot(dataEEG);
% set(gca,'XTick',1:nsec*fs:N); 
% set(gca,'XTickLabel', 0:nsec:N/fs);
% xlabel('Time (seconds)');
% ylabel('Amplitude (\muV)');
% legend('Normalized data');
% subplot(2,1,2)
% % plot(ypred_discrete);
% hold on
% plot(ypred_hybrid);
% ylim([-1 1]);
% xlim([0 N]);
% set(gca,'XTick',1:nsec*fs:N); 
% set(gca,'XTickLabel', 0:nsec:N/fs);
% xlabel('Time (seconds)');
% ylabel('Amplitude (\muV)');
% legend('Discrete Adaptive Autoregressive model','Hybrid Adaptive Autoregressive model');
% 
% second_der_hybrid=diff(ypred_hybrid,2);
% second_der_discrete=diff(ypred_discrete,2);
% % 
% % smoothness_hybrid=trapz(second_der_hybrid.^2);
% smoothness_discrete=mean(trapz(second_der_discrete.^2));
% 
% fitadaptiveNMSEhybrid=goodnessOfFit(dataEEG,ypred_hybrid,'NMSE');
% fitadaptiveNMSEdiscrete=goodnessOfFit(dataEEG,ypred_discrete,'NMSE');
% 
% %% Spectrogram estimation
psdresolution=200;
maxfreq=30;
minfreq=0.00001;
% specgram_hybrid = HybridSpectrogram(a_hybrid, minfreq, maxfreq, R,'fs', 250, 'psdres', psdresolution); 
specgram_discrete = HybridSpectrogram(a_discrete, minfreq, maxfreq, R,'fs', 250, 'psdres', psdresolution); 
% 
figure
% subplot(2,1,1)
% imagesc(specgram_hybrid.timestamp, specgram_hybrid.w, specgram_hybrid.specgram);
% xlim([0 size(dataEEG,1)])
% set(gca,'XTick',1:nsec*fs:N); 
% set(gca,'XTickLabel', 0:round(nsec/60):N/fs/60);
% xlabel('Time (Minutes)');
% caxis([-15 10]);
% colorMap=jet;
% colormap(colorMap);
% colorBar=colorbar;
% ylabel(colorBar,'Power/Frequency (dB/Hz)');
% axis xy;
% ylabel('Freq (Hz)');
% title('Hybrid Kalman Filter Spectrogram');

% subplot(2,1,2)
imagesc(specgram_discrete.timestamp, specgram_discrete.w, specgram_discrete.specgram);
xlim([0 size(dataEEG,1)])
set(gca,'XTick',1:nsec*fs:N); 
set(gca,'XTickLabel', 0:round(nsec/60):N/fs/60);
xlabel('Time (Minutes)');
caxis([-15 10]);
colorMap=jet;
colormap(colorMap);
colorBar=colorbar;
ylabel(colorBar,'Power/Frequency (dB/Hz)');
axis xy;
ylabel('Freq (Hz)');
title('Discrete Kalman Filter Spectrogram');
%% Plot spectrum every nsec seconds
%  figure
%  Legend=cell(N/indexstep,1);
%  Ffs=psdresolution/maxfreq; % sampling rate of the frequencies (index/Hz) 
% for i=1:N/indexstep
%     Legend{i}=strcat('Spectrum at t= ', num2str(i*nsec),' ',' seconds');
% end
% for i=1:N/indexstep
%     plot(specgramcont.specgram(:,indexstep*i))
%     xlabel('Frequency (Hz)');
%     ylabel('Power (dB)');
%     set(gca,'XTick',1:Ffs*(maxfreq/10):psdresolution); 
%     set(gca,'XTickLabel',1:maxfreq/10:maxfreq); % show frequency every maxfreq/10 Hz
%     hold on
% end 
% legend(Legend)
%% Plot coefficient evolution
% a_discrete=a_discrete';
% a_hybrid=a_hybrid';
% 
% 
% second_der_hybrid_coeffs=diff(a_hybrid,2);
% second_der_discrete_coeffs=diff(a_discrete,2);
% 
% smoothness_hybrid_coeffs=mean(trapz(second_der_hybrid_coeffs.^2));
% smoothness_discrete_coeffs=mean(trapz(second_der_discrete_coeffs.^2));

% figure
% for i=1:p
%     subplot(ceil(p/2),2,i)
%     plot(a_discrete(:,i))
%     xlim([0 size(dataEEG,1)])
%     set(gca,'XTick',1:nsec*fs:N); 
%     set(gca,'XTickLabel', 0:round(nsec/60):N/fs);
%     xlabel('Time (min)');
%     ylabel('Amplitude');
%     title(['a',num2str(i),' ','discrete estimation']);
% end

% figure
% for i=1:p
%     subplot(ceil(p/2),2,i)
%     plot(a_hybrid(:,i))
%     xlim([0 size(dataEEG,1)])
%     set(gca,'XTick',1:nsec*fs:N); 
%     set(gca,'XTickLabel', 0:round(nsec/60):N/fs);
%     xlabel('Time (min)');
%     ylabel('Amplitude');
%     title(['a',num2str(i),' ','hybrid estimation']);
% end