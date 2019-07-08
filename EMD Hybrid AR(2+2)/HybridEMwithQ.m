function [R_EM,Q_EM,L_EM,CL_EM] = HybridEMwithQ(y,p,P0,maxIter,Rstart,Qstart,fs,exp_smoothing,tol)
% This function applies an EM algorithm to estimate observation noise
% variance R (scalar) and process noise covariance matrix Q (pxp) of 
% hybrid state space model
%% Minibatch size to compute EM parameters
N=size(y,1);
%% Initialization of matrices to store parameter evolution
R=zeros(maxIter,1);
Q=zeros(p,p,maxIter);
L=zeros(maxIter,1); 
MarginalL=zeros(maxIter,1);
%% Calculate expected value of complete log-likelihood using initialized parameters (1st E-Step)
Q(:,:,1)=Qstart; % fix an initial process covariance matrix
R(1)=Rstart;
[a,E,Psmooth,Ppreds,H,S]=HybridKalmanEM(y,p,P0,R(1),Q(:,:,1),fs,exp_smoothing);% calculate state estimations

%%%% Complete log-likelihood %%%%
sumR1=0; % E[a(k)|z]
sumR2=0; % E[a(k)*a(k)'|z]
for k=1:N
    sumR1=sumR1+y(k)*H(k,:)*a(:,k);
    sumR2=sumR2+H(k,:)*(Psmooth(:,:,k)+a(:,k)*a(:,k)')*H(k,:)';
end
sumQ1=0; % E[a(k)*a(k)'|z]
sumQ2=0; % E[a(k)*a(k-1)'|z]
sumQ3=0; % E[a(k-1)*a(k-1)'|z]
for i=2:N
    sumQ1=sumQ1+Psmooth(:,:,i)+a(:,i)*a(:,i)'; 
    sumQ2=sumQ2+S(:,:,i-1)*Psmooth(:,:,i)+a(:,i)*a(:,i-1)'; 
    sumQ3=sumQ3+Psmooth(:,:,i-1)+a(:,i-1)*a(:,i-1)'; 
end

SSE1=sum(y.^2)-2*sumR1+sumR2;
SSE2=trace(Q(:,:,1)\(sumQ1-2*sumQ2+sumQ3));

L(1)=-N/2*log(2*pi*R(1))-((N-1)*p/2)*log(2*pi)-((N-1)/2)*log(det(Q(:,:,1)))-1/(2*R(1))*SSE1-1/2*SSE2; 
%%%% Marginal log-likelihood %%%%
sumML1=0;
sumML2=0;
for k=1:N
    temp=R(1)+H(k,:)*Ppreds(:,:,k)*H(k,:)';
    sumML1=sumML1+log(temp);
    sumML2=sumML2+(E(k)^2)/(temp);
end

MarginalL(1)=-N/2*log(2*pi)-1/2*sumML1-1/2*sumML2;
%% EM Iterations Initialization
changeML=1e+10; % initialize change in likelihood (has to be bigger than tol so the while loop starts)
j=2;
nIter=0;
%% EM Iterations
while nIter<maxIter && abs(changeML)>tol 

    %% M-Step
    %%% Observation Noise Variance Update
    R(j)=(1/N)*SSE1;
    if R(j)<0.3
        R(j)=R(1);
    end
    %%% Process Noise Covariance Update
    Q(:,:,j)=(1/(N-1))*(sumQ1-sumQ2);
     if(eig(Q(:,:,j)))<0
        warning('Negative covariance matrix eigenvalues');
    end
Q(:,:,j);
    %% E-Step
    [a,E,Psmooth,Ppreds,H,S]=HybridKalmanEM(y,p,P0,R(j),Q(:,:,j),fs,exp_smoothing); % calculate state estimations
    
    %%%% Complete log-likelihood %%%%
    sumR1=0; % E[a(k)|z]
    sumR2=0; % E[a(k)*a(k)'|z]
    for k=1:N
        sumR1=sumR1+y(k)*H(k,:)*a(:,k);
        sumR2=sumR2+H(k,:)*(Psmooth(:,:,k)+a(:,k)*a(:,k)')*H(k,:)';
    end
    sumQ1=0; % E[a(k)*a(k)'|z]
    sumQ2=0; % E[a(k)*a(k-1)'|z]
    sumQ3=0; % E[a(k-1)*a(k-1)'|z]
    for i=2:N
        sumQ1=sumQ1+Psmooth(:,:,i)+a(:,i)*a(:,i)';
        sumQ2=sumQ2+S(:,:,i-1)*Psmooth(:,:,i)+a(:,i)*a(:,i-1)'; 
        sumQ3=sumQ3+Psmooth(:,:,i-1)+a(:,i-1)*a(:,i-1)';
    end

    SSE1=sum(y.^2)-2*sumR1+sumR2;
    SSE2=trace(Q(:,:,j)\(sumQ1-2*sumQ2+sumQ3));

    L(j)=-N/2*log(2*pi*R(j))-((N-1)*p/2)*log(2*pi)-(N-1)/2*log(det(Q(:,:,j)))-1/(2*R(j))*SSE1-1/2*SSE2;

    %%%% Marginal log-likelihood %%%%
    sumML1=0;
    sumML2=0;
    for k=1:N
        temp=R(1)+H(k,:)*Ppreds(:,:,k)*H(k,:)';
        sumML1=sumML1+log(temp);
        sumML2=sumML2+(E(k)^2)/(temp);
    end

    MarginalL(j)=-N/2*log(2*pi)-1/2*sumML1-1/2*sumML2;

    %% Change in log-likelihood between the iterations
    changeML=MarginalL(j)-MarginalL(j-1);
    if changeML<0
        warning('Non monotonous evolution of marginal log-likelihood');
    end
    %% Update indices
    j=j+1;
    nIter=nIter+1;
end
%% Monitor evolution of EM algorithm iterations
%%% Observation Noise Variance
R = R(R ~= 0); % remove 0 entries in case the EM stopped earlier than maxIter
% figure
% plot(R);
% title('Observation noise variance estimates in function of iterations of EM algorithm');
% xlabel('EM Iterations');
% ylabel('Observation noise variance estimates');
%    
% %%% Complete log-likelihood convergence plot in function of EM iterations 
figure
subplot(121)
L = L(L ~= 0); % remove 0 entries in case the EM stopped earlier than maxIter
[~,idxLmax]=max(L);
plot(0:nIter,L,'-p','MarkerIndices',idxLmax,...
   'MarkerFaceColor','red','MarkerSize',12); 
xlabel('EM Iterations');
ylabel('Complete log-likelihood');
title('Complete log-likelihood');
subplot(122)
semilogy(diff(L));
xlim([0 nIter]);
grid;
ylabel('Convergence')
xlabel('EM iterations')
%%% Marginal log-likelihood and convergence plot in function of EM iterations 
figure
subplot(121)
MarginalL = MarginalL(MarginalL ~= 0); % remove 0 entries in case the EM stopped earlier than maxIter
% MarginalL=MarginalL/(N*log(2));


[~,idxLmax]=max(MarginalL);
plot(0:nIter,MarginalL,'-p','MarkerIndices',idxLmax,...
   'MarkerFaceColor','red','MarkerSize',12); 
xlabel('EM Iterations');
ylabel('Marginal log-likelihood');
title('Marginal log-likelihood');

subplot(122)
semilogy(diff(MarginalL));
xlim([0 nIter]);
grid;
ylabel('Convergence')
xlabel('EM iterations')
%% Final Estimates 
Q_EM=Q(:,:,nIter+1); % nIter+1 because the 1st entry is the L corresponding to the initialized parameters 
R_EM=R(nIter+1);
L_EM=MarginalL(nIter+1); 
CL_EM=L(nIter+1);
fprintf('EM Process Noise Covariance : ') ; Q_EM
fprintf(['EM converged in  ',num2str(nIter),' ','iterations']);
end