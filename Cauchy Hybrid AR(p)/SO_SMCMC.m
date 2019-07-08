function [A,log_likelihood] = SO_SMCMC(y,p,N,R,L,c)
% This function fits an autoregressive model of order p to EEG data using 
% Sequential mARKOMonte Carlo (SMC) filtering 
%% Input Arguments
% y: time series
% p: autoregressive model order 
% N: number of particles for the SMC 
% R: observation noise variance
% c: Coefficient of covariance matrix 
% L: steps of MCMC at each time point

T=length(y); % number of samples in the time series 
%% Initialize weights for sequential importance sampling
W=zeros(N,1);
W(:)=1/N; 
A=zeros(p,T); % coeffs x time
coeffs=aryule(y,p);
%% Initialize particles 
particles=zeros(N,p);
for j=1:size(particles,1)
    particles(j,1:p)=-coeffs(2:end);
end
%% Sequential Importance Resampling
% k=0;
% index=0;

for i=1:p
    A(:,i)=-coeffs(2:end); 
end
    log_likelihood=0;


    for i=p+1:T  % iterate accross observations 
              
        %%%%% Posterior approximation by particles and weights %%%%%
        H= y(i-1:-1:i-p)'; % Observation matrix
        particles_prev=particles;
        W_prev=W;
        particles(1:N,1:p) = cauchyrnd(particles_prev(1:N,1:p),c,N,p);
        W=W_prev.*(normpdf(y(i),H*particles',R)'); % update weights for importance sampling
        W=W./sum(W); % Normalize importance weights

       
         ESS=1/(sum(W(:).^2)); % Efficient sample size 
         if ESS<N/3
            particles=datasample(particles,N,'Weights',W(:),'Replace',true); % resample N particles with replacement
            W(:)=1/N; % Reinitialize weights
         end
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
         %%%%% Sample state from posterior density and compute logL for l=0 %%%%%
         
        particles_temp=datasample(particles,N,'Weights',W(:),'Replace',true); % Sample from posterior density
        tempA=mean(particles_temp,1);
        A(:,i)=tempA;
        W_temp=W.*(normpdf(y(i),H*particles_temp',R)');
        log_likelihood=log_likelihood+log(sum(W_temp));
        factor_temp=sum(W_temp);
        W_temp=W_temp./factor_temp;
      
        for l=1:L % MCMC steps
            particles_sampled=datasample(particles,N,'Weights',W(:),'Replace',true); % Sample from posterior density
            marginal_prev=normpdf(y(i),H*particles_temp',R)';
            W_proposed=W_temp.*(normpdf(y(i),H*particles_sampled',R)')./marginal_prev;
            proposed_log_likelihood=log_likelihood-log(sum(W_temp))+log(sum(W_proposed));% Proposed log-likelihood
            factor_proposed=sum(W_proposed);
            W_proposed=W_proposed/factor_proposed; 
            
%             u=rand();
%             u<proposed_log_likelihood/log_likelihood-1e-3
change=proposed_log_likelihood-log_likelihood;
            if   change>0

                particles_temp=particles_sampled;
%              disp('MCMC made it better!')
                W_temp=W_proposed;
                W=W_proposed;
                tempA=mean(particles_sampled,1)';
                A(:,i)=tempA;
                log_likelihood=proposed_log_likelihood;
            else
                A(:,i)=tempA; % the previous of this MCMC step for the given time
            end
        end
       
    end
end