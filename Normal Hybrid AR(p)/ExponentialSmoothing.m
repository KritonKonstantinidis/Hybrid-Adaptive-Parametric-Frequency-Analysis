function [asmooth] = ExponentialSmoothing(a)
[p,N]=size(a); % ar order
asmooth = a;
    C=1;
    c=zeros(p,1);
    for n = 1+p:1:N-1
        b=asmooth(:,n-1);
        d=a(:,n);
        asmooth(:,n) = (ones(p,1)-c).*b + c.*d;
        c=(C*(d-b).^2)./(ones(p,1)+C*(d-b).^2);
    end
end

