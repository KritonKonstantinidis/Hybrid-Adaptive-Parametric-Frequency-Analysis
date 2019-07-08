function [a,E,Psmooth,Ppreds,Hs,S] = DiscreteKalmanEM(y,p,P0,sigma2,Q,smoothing)
%% Input Arguments
% y:time series 
% p:order of the continuous autoregressive model
% F: state transition matrix 
% P0:initial state covariance matrix (px2)
% sigma2:variance of the white observation noise (scalar)
% Q:Covariance matrix of the multivariate gaussian white process noise (pxp)
% fs:sampling rate of the discrete observations 
% smoothing: boolean, smooth the estimated coefficients or not 
%% Output Arguments
%a:estimated coefficients 
%E:residuals between estimates and observations
%% Check input arguments
if smoothing ~=1 && smoothing ~=0
    error('Please input smoothing (1) or not smoothing (0)');
end
%% Initialize matrices to be used in the filtering 
N=size(y,1); % Size of the time series vector
apreds=zeros(p,N); % matrix of the predicted estimates of the coefficients
a=zeros(p,N);% matrix of the estimates of the coefficients
E= zeros(N,1);% vector with the residual error
Ppred = eye(p,p); % A priori state covariance matrix
S=zeros(p,p,N);
ypred=zeros(N,1);
Hs=zeros(N,p); % store observation matrices
Ppreds=zeros(p,p,N); % store predicted state covariance matrices at every step
Pfilt=zeros(p,p,N); % store updated state covariance matrices at every step

%% Normalize signal
% y = y./abs(hilbert(y)); % to remove amplitude bias (observation matrix
% depends on past observations, so does the Kalman gain so if they have a large value then the
% estimates will have a bigger variance) 
threshE = max(abs(y));
%% Initialize the state (coefficients) using Yule-Walker equations
ainit = aryule(y,p);
% leave out the a0 (a0=1), use minus because aryule supposes all ar terms in
% left side while in the model only y[k] is in the left side
a(:,p) = -ainit(2:end);
%% Calculate steady state for Kalman Gain and state covariance matrix
H = y(p:-1:1)';
     for i=1:15
      K = (Ppred*H')./(H*Ppred*H' + sigma2);
      Ppred=P0+Q;
      P0 = Ppred - K*H*Ppred; 
    end
%% Hybrid Kalman filtering 
    for i = 1+p:N
        
       H= y(i-1:-1:i-p)'; % Observation matrix 
       K = (Ppred*H')./(H*Ppred*H'+sigma2); % Time-varying Kalman Gain
       
      % Prediction of the state vector a(k-1|k-1)-> a(k|k-1) 
       apred=a(:,i-1);
       apreds(:,i)=apred; % store
      % Prediction of the covariance matrix P(k-1|k-1)-> P(k|k-1) 
       Ppred=P0+Q;
       S(:,:,i) = P0/(Ppred + Q);

      % Innovations
      E(i) = y(i)-H*apred;
%     It is assumed that the changes of the AAR-parameters within one 
%     interation are smaller than the estimation error.
      E(i) = min(threshE, E(i)); 
      E(i) = max(-threshE, E(i));

      % Update state vector estimate a(k-1|k)->a(k|k)
      a(:,i) = apred + K*E(i); 
     
      % Update covariance matrix  P(k-1|k)-> P(k|k) 
      P0 = Ppred - K*H*Ppred; 
      Pfilt(:,:,i)=P0; % store
      ypred(i)=H*a(:,i);
    end
    %% Kalman Smoothing
    aksmooth = zeros(size(a));
    aksmooth(:,end) = a(:,end);
    aksmooth(:,end-1) = a(:,end);
    Psmooth=zeros(p,p,N); % store smoothed state covariance matrices at every step
    Psmooth(:,:,end)=Pfilt(:,:,end); 
    Psmooth(:,:,end-1)=Pfilt(:,:,end);
    for i = N-1:-1:1+p
        aksmooth(:,i) = a(:,i) + S(:,:,i)*(aksmooth(:,i+1) - apreds(:,i+1)); % a(k|N)
        Psmooth(:,:,i) = Pfilt(:,:,i) + S(:,:,i)*(Psmooth(:,:,i+1) - Ppreds(:,:,i+1))*S(:,:,i)'; % P(k|N)
    end
%% Smoothing of the estimates
if smoothing==1
  asmooth = a;
  C=1;
  c=zeros(p,1);
  for n = 1+p:1:N-1
      b=asmooth(:,n-1);
      d=a(:,n);
      asmooth(:,n) = (ones(p,1)-c).*b + c.*d;
      c=(C*(d-b).^2)./(ones(p,1)+C*(d-b).^2);
  end
  a = asmooth;
%% Calculate the error after smoothing
  for i = 1+p:N
    E(i) = y(i) - y(i-1:-1:i-p)'*a(:,i);
  end
end
%% Fill in the inital estimates
    for k = 1:p+1
      a(:,k) = a(:,p+2);
      E(i) = E(p+2);
    end
end
