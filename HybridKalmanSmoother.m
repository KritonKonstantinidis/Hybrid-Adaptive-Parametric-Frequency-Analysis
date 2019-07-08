function [aksmooth,Psmooth,E_smooth] = HybridKalmanSmoother(y,p,afilt,apreds,Pfilt,Ppreds,S)
% This function performs the smoothing part of the Hybrid Kalman Filter
%% Input Arguments
% y: Time series 
% p: Order of the continuous autoregressive model
% afilt: Filtered coefficients from HybridKalmanFilter
% apreds: Predicted coefficients from HybridKalmanFilter
% Pfilt: Filtered state covariance matrices from HybridKalmanFilter
% Ppreds: Predicted state covariance matrices from HybridKalmanFilter
% S: Smoothing matrices from HybridKalmanFilter
%% Output Arguments
% aksmooth: Smoothed coefficients
% Psmooth: Covariance Matrices
% E_pred: Residuals between estimates and observations
%% Initialize matrices to be used in the filtering 
N=size(y,1); % Size of the time series vector
E_smooth= zeros(N,1);% Vector with the residual error calculated with smoothed coefficients
%% Kalman smoothing of the estimates
    aksmooth = zeros(size(afilt));
    aksmooth(:,end) = afilt(:,end);
    aksmooth(:,end-1) = afilt(:,end);
    Psmooth=zeros(p,p,N); % store smoothed state covariance matrices at every step
    Psmooth(:,:,end)=Pfilt(:,:,end); 
    Psmooth(:,:,end-1)=Pfilt(:,:,end);
    for i = N-1:-1:1+p
        aksmooth(:,i) = afilt(:,i) + S(:,:,i)*(aksmooth(:,i+1) - apreds(:,i+1)); % a(k|N)
        Psmooth(:,:,i) = Pfilt(:,:,i) + S(:,:,i)*(Psmooth(:,:,i+1) - Ppreds(:,:,i+1))*S(:,:,i)'; % P(k|N)
    end
%% Calculate the error after smoothing
    for i = 1+p:N
        E_smooth(i) = y(i) - y(i-1:-1:i-p)'*aksmooth(:,i);
    end
%% Fill in the inital estimates
    for k = 1:p+1
        aksmooth(:,k) = aksmooth(:,p+2);
        E_smooth(i) = E_smooth(p+2);
    end
end