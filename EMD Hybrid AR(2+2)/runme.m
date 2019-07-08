%% Load data
close all
clear all
clc

fs = 250;

data = load('SED10.mat');
data=data.data;
channel=1; % Frontal Channel
dataEEG = data(channel,600*fs+1:1880*fs);
% dataEEG=dataEEG(1:5000);
tic
EMD_Hybrid_Spectrogram(dataEEG);
toc
% N=size(dataEEG,2); % number of observations
% % 
% duration=N/fs; % duration in seconds
% nsec=duration/5; % show x axis ticks/plot instantaneous spectra every duration/5 seconds
% indexstep=nsec*fs; 
% 
% % %% EMD
% imf= emd1(dataEEG,10)';
% imf(:,1)=[];
% % % imf= eemd(dataEEG,6,10,0.1)';
% % 
% mean_frequency=meanfreq(imf,fs);
% id_frequency=mean_frequency;
% for i=1:size(id_frequency,2)
%     if mean_frequency(i)<5
%         id_frequency(i)=0;
%     elseif mean_frequency(i)<17
%                 id_frequency(i)=1;
%     else
%         id_frequency(i)=2;
%     end
% end
% % 
% slow=sum(imf(:,id_frequency==0),2);
% alpha=sum(imf(:,id_frequency==1),2);
%%
% figure
% subplot(2,1,1)
% plot(slow)
% set(gca,'XTick',1:nsec*fs:N); 
% set(gca,'XTickLabel', 0:nsec:N/fs);
% xlabel('Time (seconds)');
% ylabel('Amplitude (\muV)');
% title('Slow reconstructed oscillations');
% 
% subplot(2,1,2)
% plot(alpha)
% set(gca,'XTick',1:nsec*fs:N); 
% set(gca,'XTickLabel', 0:nsec:N/fs);
% xlabel('Time (seconds)');
% ylabel('Amplitude (\muV)');
% title('Alpha reconstructed oscillations');
% 
% D=size(imf,2); % Number of IMFs
% figure
% for i=1:D
% subplot(D,1,i)
% plot(imf(:,i))
% set(gca,'XTick',1:nsec*fs:N); 
% set(gca,'XTickLabel', 0:nsec:N/fs);
% xlabel('Time (seconds)');
% ylabel('Amplitude (\muV)');
% end
% % % 
% % % 
% reconstructed=slow+alpha;
% % % % 
% lowEnd = 0.5; 
% highEnd = 50; 
% filterOrder = 4; 
% [b, a] = butter(filterOrder, [lowEnd highEnd]/(fs/2)); % Generate filter coefficients
% dataEEG=filtfilt(b, a, dataEEG); % Apply filter to data using zero-phase filtering
% figure
% plot(dataEEG)
% hold on
% plot(reconstructed);
% set(gca,'XTick',1:nsec*fs:N); 
% set(gca,'XTickLabel', 0:nsec:N/fs);
% xlabel('Time (seconds)');
% ylabel('Amplitude (\muV)');
% legend('Original EEG data','Reconstructed');

% legend('Original EEG data','Reconstructed');
% % %% Hybrid Kalman Filtering to slow and alpha time domain signals
% % slow=slow./max(abs(slow));
% % alpha=alpha./max(abs(alpha));
% % 
% tol=1e-3;
% nIter=100;
% p=2;
% P0=eye(p);
% Rstart=0.5;
% Qstart=eye(p);
% 
% minibatch_slow=slow(1:1250);
% [Rslow,Qslow,~]=HybridEMwithQ(minibatch_slow,p,P0,nIter,Rstart,Qstart,fs,0,tol);
% minibatch_alpha=alpha(1:1250);
% [Ralpha,Qalpha,~]=HybridEMwithQ(minibatch_alpha,p,P0,nIter,Rstart,Qstart,fs,0,tol);
% % % Rslow=0.5;
% % % Qslow=1e-4*eye(p);
% % 
% a_slow=HybridKalmanFilter(slow,p,P0,Rslow,Qslow,fs,0);
% a_alpha=HybridKalmanFilter(alpha,p,P0,Ralpha,Qalpha,fs,0);
% % 
% % 
% % %% Hybrid Spectrograms 
% % % figure
% % %     
% % % specgram = HybridSpectrogram(a_slow, 0.0001, 30, 0.001,'fs', 250, 'psdres', 200); 
% % % % subplot(2,1,1)
% % % imagesc(specgram.timestamp, specgram.w, specgram.specgram);
% % % xlim([0 size(dataEEG,2)])
% % % set(gca,'XTick',1:nsec*fs:N); 
% % % set(gca,'XTickLabel', 0:round(nsec/60):N/fs/60);
% % % xlabel('Time (Minutes)');
% % % caxis([-15 10]);
% % % colorMap=jet;
% % % colormap(colorMap);
% % % colorBar=colorbar;
% % % ylabel(colorBar,'Power/Frequency (dB/Hz)');
% % % axis xy;
% % % ylabel('Freq (Hz)');
% % % title('EMD hybrid spectrogram for slow');
% % % figure
% % % specgram = HybridSpectrogram(a_alpha, 0.0001, 30, 0.001,'fs', 250, 'psdres', 200); 
% % % subplot(2,1,2)
% % % imagesc(specgram.timestamp, specgram.w, specgram.specgram);
% % % xlim([0 size(dataEEG,2)])
% % % set(gca,'XTick',1:nsec*fs:N); 
% % % set(gca,'XTickLabel', 0:round(nsec/60):N/fs/60);
% % % xlabel('Time (Minutes)');
% % % caxis([-15 10]);
% % % colorMap=jet;
% % % colormap(colorMap);
% % % colorBar=colorbar;
% % % ylabel(colorBar,'Power/Frequency (dB/Hz)');
% % % axis xy;
% % % ylabel('Freq (Hz)');
% % % title('EMD hybrid spectrogram for alpha');
% % % %
% % 
% % % figure
% % % plot(specgram.specgram(:,1))
% % %     xlabel('Frequency (Hz)');
% % %     ylabel('Power (dB)');
% % 
% % % figure
% specgram = SpectrogramEMD(a_slow,a_alpha, 0.0001, 30, 0.001,'fs', 250, 'psdres', 200); 
% imagesc(specgram.timestamp, specgram.w, specgram.specgram);
% xlim([0 size(dataEEG,2)])
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
% title('EMD hybrid spectrogram');
