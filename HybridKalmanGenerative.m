function [Z,A] = HybridKalmanGenerative(R,Q,K,p)
% This funciton uses a kalman filter with observation noise variance R and
% process noise covariance Q to generate K observations of an
% autoregressive model of order p

%% Initialize data structures 
Z=zeros(K,1);
A=zeros(p,K);
%% Sample initial p hidden states
rng default
mu_a0=zeros(1,p);
P0=1e-3*eye(p);
A(:,1:p)=mvnrnd(mu_a0,P0,p); % Draw initial p coefficient vectors
%% Generate initial p observations
for i=1:p
    Z(i)= randn(1,1);
end
%% Generate the next K-p observations 
for i = 1+p:K  
    H= Z(i-1:-1:i-p)'; % Observation matrix 
    A(:,i)=mvnrnd(A(:,i-1)',Q,1); % sample next hidden state
    Z(i)= mvnrnd(H*A(:,i),R,1);
end

end
