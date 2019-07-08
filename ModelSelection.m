function [p,Q,R,P0] = ModelSelection(pmin,pmax,criterion,y,minibatch_size,Rstart,Qstartcoef,fs,exp_smoothing,maxIter,tol_EM)
% This function calculates AIC/BIC for autoregressive model selection,
% selects the order accordingly and returns the necessary parameters for
% the hybrid Kalman filtering

AICprevious=1e+15; 
BICprevious=1e+15;
L_EM=zeros(pmax-pmin+1,1);
CL_EM=zeros(pmax-pmin+1,1);
AIC=zeros(pmax-pmin+1,1);
BIC=zeros(pmax-pmin+1,1);
j=1;
minibatch=y(round(numel(y)/4):round(numel(y)/4)+minibatch_size);
% minibatch=y(1:minibatch_size);

for i=pmin:pmax
    %%% EM for each order
    Qstart=Qstartcoef*eye(i);
    P0test=eye(i);
    
    %[R_EM,Q_EM,L_EM(j),CL_EM(j)]=HybridEMwithQ(minibatch,i,P0test,maxIter,Rstart,Qstart,fs,exp_smoothing,tol_EM);
    [R_EM,Q_EM,L_EM(j)]=HybridEM(minibatch,i,P0test,maxIter,Rstart,Qstart,fs,exp_smoothing,tol_EM);
    
    if(eig(Q_EM))<0
        warning('Negative covariance matrix eigenvalues');
    end
    %%% Information criteria to decide optimal order using results of EM
    if criterion==0
        BIC(j)=log(minibatch_size)*i-2*L_EM(j);
        if BIC(j)<BICprevious
            p=i; % chosen order
            P0=eye(p);
            Q=Q_EM;
            R=R_EM;
            BICprevious=BIC(j);
        end
    end
    if criterion==1
        AIC(j)=-2*L_EM(j)+2*i;
        if AIC(j)<AICprevious
            p=i; % chosen order
            P0=eye(p);
            Q=Q_EM;
            R=R_EM;
            AICprevious=AIC(j);
        end
    end
    j=j+1;
end
%% Information Criteria Monitoring
index=pmin:pmax;
%%% Marginal log-likelihood %%%
figure
[~,idxLmax]=max(L_EM);
plot(index,L_EM,'-p','MarkerIndices',idxLmax,...
    'MarkerFaceColor','red',...
    'MarkerSize',12);
title('Maximum marginal log likelihood in function of model order');
xlabel('Model Order');
xlim([pmin pmax]);
set(gca,'XTick',pmin:pmax); 
set(gca,'XTickLabel',pmin:pmax);
ylabel('Marginal log-likelihood');
%%% Complete log-likelihood %%%
% figure
% [~,idxCLmax]=max(CL_EM);
% plot(index,CL_EM,'-p','MarkerIndices',idxCLmax,...
%     'MarkerFaceColor','red',...
%     'MarkerSize',12);
% title('Complete log likelihood in function of model order');
% xlabel('Model Order');
% xlim([pmin pmax]);
% set(gca,'XTick',pmin:pmax); 
% set(gca,'XTickLabel',pmin:pmax);
% ylabel('Complete log-likelihood');
%%% Information Criteria %%%
if criterion==0
    figure
    [~,idxminBIC]=min(BIC);
    plot(index,BIC,'-p','MarkerIndices',idxminBIC,...
    'MarkerFaceColor','red',...
    'MarkerSize',12);
    title('BIC in function of model order');
    xlabel('Model Order');
    xlim([pmin pmax]);
    set(gca,'XTick',pmin:pmax); 
    set(gca,'XTickLabel',pmin:pmax);
    ylabel('BIC');
end

if criterion==1
    figure
    [~,idxminAIC]=min(AIC);
    plot(index,AIC,'-p','MarkerIndices',idxminAIC,...
    'MarkerFaceColor','red',...
    'MarkerSize',12);
    title('AIC in function of model order');
    xlabel('Model Order');
    xlim([pmin pmax]);
    set(gca,'XTick',pmin:pmax); 
    set(gca,'XTickLabel',pmin:pmax);
    ylabel('AIC');
end
end