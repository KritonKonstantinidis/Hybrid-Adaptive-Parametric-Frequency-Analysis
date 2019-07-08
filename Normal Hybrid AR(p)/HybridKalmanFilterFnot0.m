function [afilt,E_pred] = HybridKalmanFilterFnot0(y,p,F,P0,R,Q,fs)
%% Descpiption
% This function performs hybrid Kalman Filtering to fit an adaptive autoregressive model
% of order p to a time series 
%% Input Arguments
% y:Time series 
% p:Order of the continuous autoregressive model
% F:State transition matrix 
% P0:Initial state covariance matrix (px2)
% R:Variance of the white observation noise (scalar)
% Q:Covariance matrix of the multivariate gaussian white process noise (pxp)
% fs:Sampling rate of the discrete observations 
% exp_smoothing:Boolean, exponential smoothing of the estimated coefficients or not 
% kalman_smoothing:Boolean, kalman smoothing of the estimated coefficients or not
%% Output Arguments
% a:Estimated coefficients (smoothed or not)
% E:Residuals between estimates and observations
% P:Smoothed process covariance matrices
% Hs:Observation matrices
% Sm:Kalman smoothing matrices
% y: normalized time series
%% Initialize matrices to be used in the filtering 
N=size(y,1); % Size of the time series vector
apreds=zeros(p,N);
afilt=zeros(p,N);% matrix of the updated estimates of the coefficients
E_pred= zeros(N,1);% vector with the residual error calculated with predicted coefficients
Ppred = eye(p,p); % A priori state covariance matrix
Ppreds=zeros(p,p,N);
Pfilt=zeros(p,p,N); % stored updated state covariance matrices at every step
%% Initialize the state (coefficients) using Yule-Walker equations
aInit = aryule(y,p);
% leave out the a0 (a0=1), use minus because aryule supposes all ar terms in
% left side while in the model only y[k] is in the left side
afilt(:,p) = -aInit(2:end);
%% Calculate steady state for Kalman gain and state covariance matrix
H = y(p:-1:1)';
for i=1:15
    K = (Ppred*H')./(H*Ppred*H' + R);
    [~,Ppred] = ode45(@(t,P)mLyapunov(t, P, F, Q), [0 1/fs], P0);% Prediction of the covariance matrix P(k-1|k-1)-> P(k|k-1)
    Ppred=Ppred(end,:); % the prediction at the last step of the integration
    Ppred=reshape(Ppred,[p,p]); % reshape the covariance matrix into a pxp matrix
    P0 = Ppred - K*H*Ppred; 
end
%% Hybrid Kalman filtering
for i = 1+p:N
        
    H= y(i-1:-1:i-p)'; % Observation matrix 
    K = (Ppred*H')./(H*Ppred*H'+R); % Time-varying Kalman Gain (scalar division) for each state variable

    % Prediction of the state vector a(k-1|k-1)-> a(k|k-1) 
    [~,apred] = ode45(@(t,a)mDifferential(t, a, F), [0 1/fs], afilt(:,i-1));
    apred=apred(end,:)'; % the prediction at the last time step of the integration
    apreds(:,i)=apred; % store
    
    % Prediction of the covariance matrix P(k-1|k-1)-> P(k|k-1) 
    [~,Ppred] = ode45(@(t,P)mLyapunov(t, P, F, Q), [0 1/fs], P0);
    Ppred=Ppred(end,:); % the prediction at the last step of the integration
    Ppred=reshape(Ppred,[p,p]); % reshape the covariance matrix into a pxp matrix
    Ppreds(:,:,i)=Ppred; % store
    
    E_pred(i) = y(i)-H*apred; % Innovations

    % Update state vector estimate a(k-1|k)->a(k|k)
    afilt(:,i) = apred + K*E_pred(i); 
    
    % Update covariance matrix  P(k-1|k)-> P(k|k) 
    P0 = Ppred - K*H*Ppred; 
    Pfilt(:,:,i)=P0;
end
%% Fill in the inital estimates
    for k = 1:p+1
        afilt(:,k) = afilt(:,p+2);
        E_pred(i) = E_pred(p+2);
    end
end

function dXdt = mLyapunov(~, X, A, Q)
    X = reshape(X, size(A)); %Convert from "n^2"-by-1 to "n"-by-"n"
    dXdt = A*X + X*A' +Q; %Determine derivative
    dXdt = dXdt(:); %Convert from "n"-by-"n" to "n^2"-by-1
end

function dXdt = mDifferential(~, X, A)
    dXdt = A*X; %Determine derivative
end