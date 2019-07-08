function [aksmooth,E_pred,Psmooth,Ppreds,Hs,S] = HybridKalmanEM(y,p,P0,R,Q,fs,exp_smoothing)
%% Descpiption
% This function performs hybrid Kalman Filtering and Smoothing to fit an adaptive autoregressive model
% of order p to a time series y
%% Input Arguments
% y: Time series 
% p: Order of the continuous autoregressive model
% P0: Initial state covariance matrix (px2)
% R: Variance of the white observation noise (scalar)
% Q: Covariance matrix of the multivariate gaussian white process noise (pxp)
% fs: Sampling rate of the discrete observations 
% exp_smoothing: Boolean, exponential smoothing of the estimated coefficients or not 
%% Output Arguments
% aksmooth:Estimated smoothed coefficients
% E_pred:Residuals between estimates and observations
% Psmooth:Smoothed process covariance matrices
% Hs:Observation matrices
% S:Kalman smoothing matrices
%% Check input arguments
if exp_smoothing ~=1 && exp_smoothing ~=0
    error('Please input exponential smoothing (1) or not exponential smoothing (0)');
end
%% Initialize matrices to be used in the filtering 
N=size(y,1); % Size of the time series vector
apreds=zeros(p,N); % matrix of the predicted estimates of the coefficients
afilt=zeros(p,N); % matrix of the updated estimates of the coefficients
E_pred= zeros(N,1);% vector with the residual error calculated with predicted coefficients
Ppred = eye(p,p); % A priori state covariance matrix
S=zeros(p,p,N); % kalman smoothing matrices
Ppreds=zeros(p,p,N); % store predicted state covariance matrices at every step
Pfilt=zeros(p,p,N); % store updated state covariance matrices at every step
Hs=zeros(N,p); % store observation matrices
threshE = max(abs(y));
%% Initialize the state (coefficients) using Yule-Walker equations
aInit = aryule(y,p);
% leave out the a0 (a0=1), use minus because aryule supposes all ar terms in
% left side while in the model only y[k] is in the left side
afilt(:,p) = -aInit(2:end);
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
    Hs(i,:)=H; % store
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
    
    E_pred(i) = y(i)-H*apred; % Innovations

    % It is assumed that the changes of the AAR-parameters within one 
    % iteration are smaller than the estimation error. Bound error to
    % prevent instabilities 
    E_pred(i) = min(threshE, E_pred(i)); 
    E_pred(i) = max(-threshE, E_pred(i));
      
    % Update state vector estimate a(k-1|k)->a(k|k)
    afilt(:,i) = apred + K*E_pred(i); 
    
    % Update covariance matrix  P(k-1|k)-> P(k|k) 
    P0 = Ppred - K*H*Ppred; 
    Pfilt(:,:,i)=P0; % store
    
    % Kalman smoothing matrix using the identity as transition, bar-shalom, 8.6.2-17
    S(:,:,i) = P0/(P0 + Q); 
end
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
%% Exponential smoothing of the estimates
if exp_smoothing==1
   aksmooth=ExponentialSmoothing(aksmooth);
end
%% Recalculate the error after smoothing
    for i = 1+p:N
        E_pred(i) = y(i) - y(i-1:-1:i-p)'*aksmooth(:,i);
    end
%% Fill in the inital estimates
    for k = 1:p+1
        aksmooth(:,k) = aksmooth(:,p+2);
        E_pred(i) = E_pred(p+2);
        Hs(k,:)=Hs(p+2,:);
    end
end

function dXdt = mLyapunov(~, ~, Q)
    dXdt = Q; %Determine derivative
    dXdt = dXdt(:); %Convert from "n"-by-"n" to "n^2"-by-1
end

function dXdt = mDifferential(~, ~,p)
    dXdt = zeros(p,1); %Determine derivative
end