function [afilt,apreds,Pfilt,Ppreds,S,E,ypred] = HybridKalmanFilterOnline(y,p,P0,R,Q,fs,exp_smoothing,aprev)
% This function performs the filtering part of the hybrid Kalman filter
%% Input Arguments
% y: Time series 
% p: Order of the continuous autoregressive model
% P0: Initial state covariance matrix (px2)
% R: Variance of the white observation noise (scalar)
% Q: Covariance matrix of the multivariate gaussian white process noise (pxp)
% fs: Sampling rate of the discrete observations 
% exp_smoothing: Boolean, exponential smoothing of the estimated coefficients or not 
%% Output Arguments
% afilt: Filtered coefficients
% apreds: Predicted coefficients
% Pfilt: Filtered state covariance matrices
% Ppreds: Predicted state covariance matrices
% S: Kalman smoothing matrices
%% Check input arguments
if exp_smoothing ~=1 && exp_smoothing ~=0
    error('Please input exponential smoothing (1) or not exponential smoothing (0)');
end
%% Initialize matrices to be used in the filtering 
N=size(y,1); % Size of the time series vector
apreds=zeros(p,N); % matrix of the predicted estimates of the coefficients
afilt=zeros(p,N); % matrix of the updated estimates of the coefficients
E= zeros(N,1);% vector with the residual error calculated with predicted coefficients
Ppred = eye(p,p); % A priori state covariance matrix
Ppreds=zeros(p,p,N); % store predicted state covariance matrices at every step
Pfilt=zeros(p,p,N); % store updated state covariance matrices at every step
S=zeros(p,p,N); % kalman smoothing matrices
threshE = max(abs(y));
ypred=zeros(N,1);
%% Initialize the state (coefficients) using the estimates from the previous filtering procedure
afilt(:,p) = aprev;
%% Calculate steady state for Kalman gain and state covariance matrix
H = y(p:-1:1)';
for i=1:15
    K = (Ppred*H')./(H*Ppred*H' + R);
    [~,Ppred] = ode45(@(t,P)mLyapunov(t, P, Q), [0 1/fs], P0);% Prediction of the covariance matrix P(k-1|k-1)-> P(k|k-1)
    Ppred=Ppred(end,:); % the prediction at the last step of the integration
    Ppred=reshape(Ppred,[p,p]); % reshape the covariance matrix into a pxp matrix
    P0 = Ppred - K*H*Ppred; 
end
%% Hybrid Kalman filtering
for i = 1+p:N
        
    H= y(i-1:-1:i-p)'; % Observation matrix 
    K = (Ppred*H')./(H*Ppred*H'+R); % Time-varying Kalman Gain (scalar division) for each state variable

    % Prediction of the state vector a(k-1|k-1)-> a(k|k-1) 
    [~,apred] = ode45(@(t,a)mDifferential(t, a, p), [0 1/fs], afilt(:,i-1));
    apred=apred(end,:)'; % the prediction at the last time step of the integration
    apreds(:,i)=apred; % store

    % Prediction of the covariance matrix P(k-1|k-1)-> P(k|k-1) 
    [~,Ppred] = ode45(@(t,P)mLyapunov(t, P, Q), [0 1/fs], P0);
    Ppred=Ppred(end,:); % the prediction at the last step of the integration
    Ppred=reshape(Ppred,[p,p]); % reshape the covariance matrix into a pxp matrix
    Ppreds(:,:,i)=Ppred; % store

    E_pred = y(i)-H*apred; % Innovations

    % It is assumed that the changes of the AAR-parameters within one 
    % iteration are smaller than the estimation error. Bound error to
    % prevent instabilities 
    E_pred = min(threshE, E_pred); 
    E_pred = max(-threshE, E_pred);
    
    % Update state vector estimate a(k-1|k)->a(k|k)
    afilt(:,i) = apred + K*E_pred; 
    E(i) = y(i)-H*afilt(:,i); 

    % Update covariance matrix  P(k-1|k)-> P(k|k) 
    P0 = Ppred - K*H*Ppred; 
    Pfilt(:,:,i)=P0; % store
    
    % Kalman smoothing matrix using the identity as transition
    S(:,:,i) = P0/(P0 + Q); % bar-shalom, 8.6.2-17
    
    ypred(i)=H*afilt(:,i);

end
%% Exponential smoothing of the estimates
if exp_smoothing==1
    afilt=ExponentialSmoothing(afilt);
end

%% Fill in the inital estimates
    for k = 1:p+1
        afilt(:,k) = afilt(:,p);
        E(i) = E(p);
    end
end

function dXdt = mLyapunov(~, ~, Q)
    dXdt = Q; %Determine derivative
    dXdt = dXdt(:); %Convert from "n"-by-"n" to "n^2"-by-1
end

function dXdt = mDifferential(~, ~,p)
    dXdt = zeros(p,1); %Determine derivative
end
