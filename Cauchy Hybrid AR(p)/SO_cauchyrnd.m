function [numbers] = SO_cauchyrnd(mu,c,N,p)

temp=tan(pi*(rand(N,p)-0.5))'; 
numbers= (diag(c) * temp)'+mu;

%%% Not vectorized 
% numbers=zeros(N,p);
% temp=tan(pi*(rand(N,p)-0.5))'; 
% 
% for j = 1 : N
%   numbers(j,:) = c .* temp( :, j )+mu(j,:)';
% end

end
