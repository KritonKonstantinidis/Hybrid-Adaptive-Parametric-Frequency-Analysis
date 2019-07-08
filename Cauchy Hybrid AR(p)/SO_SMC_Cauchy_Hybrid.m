function [A,log_likelihood,c] = SO_SMC_Cauchy_Hybrid(y,p,N,R,delta,a,b,R_e)
% This function fits an autoregressive model of order p to EEG data using 
% Sequential Monte Carlo (SMC) filtering 
%% Input Arguments
% y: time series
% p: autoregressive model order 
% N: number of particles for the SMC 
% R: observation noise variance
% c: Coefficient of covariance matrix 
% R_e: variance of the random walk of the variances of the cauchy distributions

T=length(y); % number of samples in the time series 

%% Initialize weights for sequential importance sampling
W=zeros(N,1);
W(:)=1/N; 
A=zeros(p,T); % coeffs x time
c=zeros(p,T);
H=zeros(1,2*p);
coeffs=aryule(y,p);
%% Initialize particles 
particles=zeros(N,2*p);
for j=1:size(particles,1)
    particles(j,1:p)=-coeffs(2:end);
    particles(j,p+1:2*p)=log(a + (b-a)*rand(1,p)); % initialize particles uniformly in [a,b] 
end
%% Sequential Importance Resampling
k=0;
% index=0;
log_likelihood=0;

for i=1:p
    A(:,i)=-coeffs(2:end); 
    c(:,i)=a + (b-a)*rand(1,p);
end
  
    for i=p+1:T         
        
         W_prev=W;
        % Observation matrix 
        H(1:p)= y(i-1:-1:i-p)'; 
        H(p+1:2*p)=zeros(1,p);
        
        for j=1:delta
            particles_prev = particles;
            particles(1:N,1:p) = SO_cauchyrnd(particles_prev(1:N,1:p),c(:,i-1)/delta,N,p);
            particles(1:N,p+1:2*p)=normrnd(particles_prev(1:N,p+1:2*p),R_e/delta,N,p);
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
            particles=datasample(particles,N,1,'Weights',W(:),'Replace',true); % resample N particles with replacement
            W(:)=1/N; % reinitialize weights
        end
        
        A(:,i)=(particles(1:N,1:p))'*W(:); % estimated coefficient is the weighted mean of the filtered density
        c(:,i)=exp((particles(1:N,p+1:2*p)'))*W(:); % update covariances of cauchy 
    end 
%k
end