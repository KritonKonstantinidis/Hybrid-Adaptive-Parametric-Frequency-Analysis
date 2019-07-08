function [numbers] = cauchyrnd(mu,c,N,p)

    numbers = c.*tan(pi*(rand(N,p)-0.5))+mu;
    
end
