function [A,log_likelihood,e,rec] = SMC_Cauchy(y,p,N,R,c)
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
e=zeros(N,1);
rec=zeros(N,1);
W(:)=1/N; % 
A=zeros(p,T); % coeffs x time
coeffs=aryule(y,p);
%% Initialize particles 
particles=zeros(N,p);
for j=1:size(particles,1)
    particles(j,:)=-coeffs(2:end);
end
%% Sequential Importance Resampling
log_likelihood=0;

for i=1:p
    A(:,i)=-coeffs(2:end); 
end
  
    for i=p+1:T       
  
        particles_prev=particles;
        W_prev=W;
        H=y(i-1:-1:i-p)'; % Observation matrix 
        
        particles = cauchyrnd(particles_prev,c,N,p);
%         particles = cauchymvrnd(particles_prev,c,p);
        W=W_prev.*(normpdf(y(i),H*particles',R)'); % update weights for importance sampling
        log_likelihood=log_likelihood+log(sum(W));

        W=W/sum(W); % Normalize importance weights
        % Resample with replacement to avoid particle degeneracy in case of
        % ESS less than threshold, use threshold to avoid particle
        % impoverishment 
        ESS=1/(sum(W(:).^2)); % Efficient sample size 
        if ESS<N/3
            particles=datasample(particles,N,'Weights',W(:),'Replace',true); % resample N particles with replacement
            W(:)=1/N; % reinitialize weights
        end
        A(:,i)=particles'*W(:); % estimated coefficient is the weighted mean of the filtered density    
        rec(i)=H*A(:,i);
        e(i)=y(i)-rec(i);

    end 
end