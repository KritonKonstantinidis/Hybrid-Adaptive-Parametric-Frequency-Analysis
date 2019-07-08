function specgram = HybridSpectrogram(A,minfreq,maxfreq,noisevar,varargin)

%% Input Arguments
% A: matrix of AR coefficients (pxn)
% maxfreq: maximum frequency until which power is calculated
% fs: sampling rate 
% index: vector of discrete time indices to do frequency decompositon
% psdres: frequency resolution
%% Output Arguments
% sgram: structure containing information to produce spectrogram
%% Frequency decomposition
  P = size(A,1); % number of coefficients
  N = size(A,2); % number of observations 

  %%% Default values 
  opt.fs = 250;
  opt.index = 1:N;
  opt.psdres = 200;
  opt.timestamp = [];
  
  opt = parsevarargin(varargin, opt); % Parse input arguments 

  if maxfreq>= opt.fs
      error('Maximum frequency must be less than the half of the sampling rate');
  end
  if isempty(opt.timestamp)
    opt.timestamp = 1:N;
  end

  N = length(opt.index);
  
  % Calculate spectrogram for frequencies below fs/2, w=2πf so f=w/(2π)
  % 0<=w<=pi 0<=w*fs/2pi<=pi*fs/2pi 0<=ffs<=fs/2 
  % we are interested in the band 0-maxfreq Hz so 0->maxfreq*pi/(fs/2)
  
  omega = linspace(minfreq*pi/(opt.fs/2), maxfreq*pi/(opt.fs/2), opt.psdres); %normalized frequencies
  
  specgram.specgram = zeros(length(omega),N); % num_freqs x observations
  specgram.w = omega.*(opt.fs./2./pi);  % actual frequencies
  specgram.timestamp = opt.timestamp;
  specgram.fs = opt.fs;
  
  for n=1:N   % for each sample
%     index = opt.index(n);
    cp = [1, -A(:,n)']'; % store the p previous coefficients, a0=1 
    
    for wi = 1:length(omega) % for each normalized frequency 
      w = omega(wi); %%% set of normalized frequencies to calculate spectrum
      z = exp(1i*w); % z transform
      powers = z.^(0:P)';
      %%%% for propofol use around 0.0015, for ketamine around 0.15
      specgram.specgram(wi,n) = noisevar/((abs(sum(powers.*cp)))^2);
    end
  end
  specgram.specgram=20*log10(specgram.specgram); % transform to dB scale
end