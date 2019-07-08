function [A,log_likelihood] = SMC_Cauchy_Hybrid(y,p,N,R,c,delta)
% This function fits an autoregressive model of order p to EEG data using 
% Sequential Monte Carlo (SMC) filtering 
%% Input Arguments
% y: time series
% p: autoregressive model order 
% N: number of particles for the SMC 
% R: observation noise variance
% c: Coefficient of covariance matrix 

T=length(y); % number of samples in the time series 

%% Initialize weights for sequential importance sampling
W=zeros(N,1);
W(:)=1/N; % 
A=zeros(p,T); % coeffs x time
coeffs=aryule(y,p);
%% Initialize particles 
particles=zeros(N,p);
for j=1:size(particles,1)
    particles(j,:)=-coeffs(2:end);
end
%% Sequential Importance Resampling
% k=0;
% index=0;
log_likelihood=0;

for i=1:p
    A(:,i)=-coeffs(2:end); 
end
  
    for i=p+1:T         
        W_prev=W;
        H= y(i-1:-1:i-p)'; % Observation matrix 
        
        for j=1:delta
            particles_prev = particles;
            particles=cauchyrnd(particles_prev,c/delta,N,p);
        end
        
        W=W_prev.*(normpdf(y(i),H*particles',R)'); % update weights for importance sampling
        log_likelihood=log_likelihood+log(sum(W));
        W=W/sum(W); % Normalize importance weights
        % Resample with replacement to avoid particle degeneracy in case of
        % ESS less than threshold, use threshold to avoid particle
        % impoverishment 
        ESS=1/(sum(W(:).^2)); % Efficient sample size 
        if ESS<N/3
%             k=k+1;
            particles=datasample(particles,N,'Weights',W(:),'Replace',true); % resample N particles with replacement
            W(:)=1/N; % reinitialize weights
        end
        A(:,i)=particles'*W(:); % estimated coefficient is the weighted mean of the filtered density        
    end 
% k
end