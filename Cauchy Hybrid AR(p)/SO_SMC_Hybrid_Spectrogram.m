function [] = SO_SMC_Hybrid_Spectrogram(y,p,N,R,B,delta,a,d)
% This function fits an autoregressive model of order p to EEG data using 
% Sequential Monte Carlo (SMC) filtering 
%% Input Arguments
% y: time series
% p: autoregressive model order 
% N: number of particles for the SMC 
% R: observation noise variance
% c: Scale parameter of the Cauchy distribution of the AR parameters
% B: Number of Bootstrap Samples
% delta: Number of discrete steps between observations 
% a,b: Interval for uniform initialization of the scale parameters

T=size(y,1); % number of observations
fs=250; % sampling rate
duration=T/fs; % duration in seconds
nsec=duration/5; % show x axis ticks/plot instantaneous spectra every duration/5 seconds

y=y./(max(abs(y))); % data normalization 

a_cauchy_hybrid=zeros(p,T,B);
likelihood_cauchy_hybrid=zeros(B,1);
for b=1:B

      [a_cauchy_hybrid(:,:,b),likelihood_cauchy_hybrid(b),c]=SO_SMC_Cauchy_Hybrid(y,p,N,R,delta,a,d);
      temp_cauchy_hybrid=a_cauchy_hybrid(:,:,b);

    sumnan_hybrid=sum(isnan(temp_cauchy_hybrid(:)));
    if(sumnan_hybrid>=1)
        disp('oops')
        a_cauchy_hybrid(:,:,b)=[];
        likelihood_cauchy_hybrid(b)=[];
    end
    
a_cauchy_hybrid=mean(a_cauchy_hybrid,3);
likelihood_cauchy_hybrid=mean(likelihood_cauchy_hybrid);

%%%% Plot the scale parameters' temporal evolution (optional)
figure
for i=1:p
    subplot(ceil(p/2),2,i)
    plot(c(i,:))
    xlim([0 size(y,1)])
    set(gca,'XTick',1:nsec*fs:N); 
    set(gca,'XTickLabel', 0:round(nsec/60):N/fs);
    xlabel('Time (min)');
    ylabel('Amplitutde');
    title(['Scale parameter for a',num2str(i)]);
end

psdresolution=200;
maxfreq=30;
minfreq=0.00001;
specgram_cauchy_hybrid = HybridSpectrogram(a_cauchy_hybrid, minfreq, maxfreq, 0.0015,'fs', 250, 'psdres', psdresolution);

figure
imagesc(specgram_cauchy_hybrid.timestamp, specgram_cauchy_hybrid.w, specgram_cauchy_hybrid.specgram);
xlim([0 size(y,1)])
set(gca,'XTick',1:nsec*fs:T); 
set(gca,'XTickLabel', 0:round(nsec/60):T/fs/60);
xlabel('Time (Minutes)');
caxis([-15 10]);
colorMap=jet;
colormap(colorMap);
colorBar=colorbar;
ylabel(colorBar,'Power/Frequency (dB/Hz)');
axis xy;
ylabel('Freq (Hz)');
title('Hybrid Sequential Monte Carlo Spectrogram');
end