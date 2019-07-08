function [] = EMD_Hybrid_Spectrogram(dataEEG)

% This function generates a time-varying spectrogram of an EEG under
% propofol using EMD and hybrid autoregressive models of order 2

%% Data Metrics
fs = 250;
N=size(dataEEG,2); % number of observations
duration=N/fs; % duration in seconds
nsec=duration/5; % show x axis ticks/plot instantaneous spectra every duration/5 seconds
%% Perform EMD
imf= emd1(dataEEG,10)';
%% Separate slow and alpha frequency bands in time domain 
mean_frequency=meanfreq(imf,fs);
id_frequency=mean_frequency;
for i=1:size(id_frequency,2)
    if mean_frequency(i)<5
        id_frequency(i)=0;
    elseif mean_frequency(i)<17
                id_frequency(i)=1;
    else
        id_frequency(i)=2;
    end
end

slow=sum(imf(:,id_frequency==0),2);
alpha=sum(imf(:,id_frequency==1),2);
%% Plot the reconstructed time-domain signals
figure

subplot(2,1,1)
plot(slow)
set(gca,'XTick',1:nsec*fs:N); 
set(gca,'XTickLabel', 0:nsec:N/fs);
xlabel('Time (seconds)');
ylabel('Amplitude (\muV)');
title('Slow reconstructed oscillations');

subplot(2,1,2)
plot(alpha)
set(gca,'XTick',1:nsec*fs:N); 
set(gca,'XTickLabel', 0:nsec:N/fs);
xlabel('Time (seconds)');
ylabel('Amplitude (\muV)');
title('Alpha reconstructed oscillations');
%% Hybrid Kalman Filtering to slow and alpha time domain signals

slow=slow./max(abs(slow));
alpha=alpha./max(abs(alpha));

tol=1e-3;
nIter=100;
p=2;
P0=eye(p);
Rstart=0.5;
Qstart=eye(p);

minibatch_slow=slow(1:1250);
[Rslow,Qslow,~]=HybridEMwithQ(minibatch_slow,p,P0,nIter,Rstart,Qstart,fs,0,tol);
minibatch_alpha=alpha(1:1250);
[Ralpha,Qalpha,~]=HybridEMwithQ(minibatch_alpha,p,P0,nIter,Rstart,Qstart,fs,0,tol);

a_slow=HybridKalmanFilter(slow,p,P0,Rslow,Qslow,fs,0);
a_alpha=HybridKalmanFilter(alpha,p,P0,Ralpha,Qalpha,fs,0);

%% Generate Hybrid Spectrogram
figure
specgram = SpectrogramEMD(a_slow,a_alpha, 0.0001, 30,'fs', 250, 'psdres', 200); 
imagesc(specgram.timestamp, specgram.w, specgram.specgram);
xlim([0 size(dataEEG,2)])
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
title('EMD hybrid spectrogram');

end