clear all
close all
clc
%%
fs=250;
%%% KETAMINE DATA
% dataEEG=load('SED974_lite.mat');
% dataEEG=dataEEG.data;

data = load('SED10.mat');
data=data.data;
channel=1; % Frontal Channel
dataEEG = data(channel,600*fs+1:1880*fs)'; % alpha oscillations
R=0.1;
dataEEG=dataEEG./(max(abs(dataEEG))); % data normalization 

N=size(dataEEG,1); % number of observations
duration=N/fs; % duration in seconds
nsec=duration/5; % show x axis ticks/plot instantaneous spectra every duration/5 seconds
indexstep=nsec*fs; 

N_particles_Cauchy=500;
c=1e-5;
B=2;

pmin=14;
pmax=14;
j=1;

for p=pmin:pmax
    a_cauchy_hybrid=zeros(p,N,B);
    likelihood_cauchy_hybrid=zeros(B,1);
    for b=1:B
        disp(['Running bootstrap sample number',' ',num2str(b)])
        [a_cauchy_hybrid(:,:,b),likelihood_cauchy_hybrid(b),e,rec]=SMC_Cauchy(dataEEG,p,N_particles_Cauchy,R,c);

%         [a_cauchy_hybrid(:,:,b),likelihood_cauchy_hybrid(b)]=SO_SMCMC(dataEEG,p,N_particles_Cauchy,R,500,c);
%         [a_cauchy_hybrid(:,:,b),likelihood_cauchy_hybrid(b),~]=SO_SMCMC_Cauchy(dataEEG,p,N_particles_Cauchy,R,5,1e-5,1e-5);
%         [a_cauchy_hybrid(:,:,b),likelihood_cauchy_hybrid(b),c]=SO_SMC_Cauchy(dataEEG,p,N_particles_Cauchy,R,1e-4,1e-6);
%        [a_cauchy_hybrid(:,:,b),likelihood_cauchy_hybrid(b),c]=SO_SMC_Cauchy_Hybrid(dataEEG,p,N_particles_Cauchy,R,delta,1e-5,1e-6,1e-3);

        temp_cauchy_hybrid=a_cauchy_hybrid(:,:,b);
    
         sumnan_hybrid=sum(isnan(temp_cauchy_hybrid(:)));
        if(sumnan_hybrid>=1)
            disp('oops')
            a_cauchy_hybrid(:,:,b)=[];
            likelihood_cauchy_hybrid(b)=[];
        end
    end
a_cauchy_hybrid=mean(a_cauchy_hybrid,3);

% likelihood_cauchy=mean(likelihood_cauchy_hybrid);
% AIC(j)=2*p-2*mean(likelihood_cauchy_hybrid);
% j=j+1;
end


%% Spectrogram 


psdresolution=200;
maxfreq=30;
minfreq=0.00001;

specgram_cauchy_hybrid = HybridSpectrogram(a_cauchy_hybrid, minfreq, maxfreq, 0.0015,'fs', 250, 'psdres', psdresolution); 
figure
imagesc(specgram_cauchy_hybrid.timestamp, specgram_cauchy_hybrid.w, specgram_cauchy_hybrid.specgram);
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
