function [R_EM,Q_EM,L_EM] = DiscreteEM(y,p,P0,maxIter,Rstart,Qstart,exp_smoothing,tol)
% This function applies an EM algorithm to estimate observation noise
% variance R (scalar) and process noise covariance matrix Q (pxp) of 
% hybrid state space model
%% Minibatch size to compute EM parameters
N=size(y,1);
%% Initialization of matrices to store parameter evolution
R=zeros(maxIter,1);
Q=zeros(p,p,maxIter);
L=zeros(maxIter,1);
%% Calculate marginal log-likelihood using initialized parameters (1st E-Step)
Q(:,:,1)=Qstart; % fix an initial process covariance matrix
R(1)=Rstart;
[a,E,Psmooth,Ppreds,H,S]=DiscreteKalmanEM(y,p,P0,R(1),Q(:,:,1),exp_smoothing);% calculate state estimations
sumR1=0; % E[a(k)|z]
sumR2=0; % E[a(k)*a(k)'|z]

for k=1:N
    sumR1=sumR1+y(k)*H(k,:)*a(:,k);
    sumR2=sumR2+H(k,:)*(Psmooth(:,:,k)+a(:,k)*a(:,k)')*H(k,:)';
end
sumQ1=0; % E[a(k)*a(k)'|z]
sumQ2=0; % E[a(k)*a(k-1)'|z]
for i=2:N
    sumQ1=sumQ1+Psmooth(:,:,i)+a(:,i)*a(:,i)'; 
    sumQ2=sumQ2+S(:,:,i-1)*Psmooth(:,:,i)+a(:,i)*a(:,i-1)'; 
end

SSE=sum(y.^2)-2*sumR1+sumR2;

%%%% Marginal log-likelihood %%%%
sumL1=0;
sumL2=0;
for k=1:N
    temp=R(1)+H(k,:)*Ppreds(:,:,k)*H(k,:)';
    sumL1=sumL1+log(temp);
    sumL2=sumL2+(E(k)^2)/(temp);
end

% % L(1)=-N/2*log(2*pi)-1/2*sum(log(R(1)+H(k,:).*Ppreds(:,:,k).*H(k,:)')+(E(k).^2)/(R(1)+H(k,:).*Ppreds(:,:,k).*H(k,:)'))

L(1)=-N/2*log(2*pi)-1/2*sumL1-1/2*sumL2;

%% EM Iterations Initialization
changeL=1e+10; % initialize change in likelihood (has to be bigger than tol so the while loop starts)
j=2;
nIter=0;
%% EM Iterations
while nIter<maxIter && abs(changeL)>tol 
    
    %% M-Step
    %%% Observation Noise Variance Update
    R(j)=(1/N)*SSE;
    if R(j)<0.3
        R(j)=R(1);
    end
    %%% Process Noise Covariance Update
    Q(:,:,j)=(1/(N-1))*(sumQ1-sumQ2);
    if(eig(Q(:,:,j)))<0
        warning('Negative covariance matrix eigenvalues');
    end
% Q(:,:,j)
    %% E-Step
    [a,E,Psmooth,Ppreds,H,S]=DiscreteKalmanEM(y,p,P0,R(j),Q(:,:,j),exp_smoothing); % calculate state estimations
    
    sumR1=0; % E[a(k)|z]
    sumR2=0; % E[a(k)*a(k)'|z]
    for k=1:N
        sumR1=sumR1+y(k)*H(k,:)*a(:,k);
        sumR2=sumR2+H(k,:)*(Psmooth(:,:,k)+a(:,k)*a(:,k)')*H(k,:)';
    end
    sumQ1=0; % E[a(k)*a(k)'|z]
    sumQ2=0; % E[a(k)*a(k-1)'|z]
    for i=2:N
        sumQ1=sumQ1+Psmooth(:,:,i)+a(:,i)*a(:,i)';
        sumQ2=sumQ2+S(:,:,i-1)*Psmooth(:,:,i)+a(:,i)*a(:,i-1)'; 
    end

    SSE=sum(y.^2)-2*sumR1+sumR2;
    
    %%%% Marginal log-likelihood %%%%
    sumL1=0;
    sumL2=0;
    for k=1:N
        temp=R(j)+H(k,:)*Ppreds(:,:,k)*H(k,:)';
        sumL1=sumL1+log(temp);
        sumL2=sumL2+(E(k)^2)/(temp);
    end

    L(j)=-N/2*log(2*pi)-1/2*sumL1-1/2*sumL2;
    %% Change in log-likelihood between the iterations
    changeL=L(j)-L(j-1);
    if changeL<0
        warning('Non monotonous evolution of marginal log-likelihood');
    end
    %% Update indices
    j=j+1;
    nIter=nIter+1;
end
%% Monitor evolution of EM algorithm iterations  
% figure
% subplot(121)
% L = L(L ~= 0); % remove 0 entries in case the EM stopped earlier than maxIter
% [~,idxLmax]=max(L);
% plot(0:nIter,L,'-p','MarkerIndices',idxLmax,...
%    'MarkerFaceColor','red','MarkerSize',12); 
% xlabel('EM Iterations');
% ylabel('Marginal log-likelihood');
% subplot(122)
% semilogy(diff(L));
% xlim([0 nIter]);
% grid;
% ylabel('Convergence')
% xlabel('EM iterations')
%% Final Estimates 
Q_EM=Q(:,:,nIter+1); % nIter+1 because the 1st entry is the L corresponding to the initialized parameters 
R_EM=R(nIter+1);
L_EM=L(nIter+1); 
% fprintf('EM Process Noise Covariance : ') ; Q_EM
% fprintf(['EM converged in  ',num2str(nIter),' ','iterations']);
end