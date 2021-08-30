function [p] = boundary_hnr(A,b,x0,N,W)
    
    d = size(A,2);
    m = size(A,1);
    x = x0;
    p=zeros(d, N);
    for i=1:N
        for j=1:W
            v = randn(d,1);
            v= v/norm(v);
            
          
            lambdas = (b - A*x) ./ (A*v);
            lambdas = 1./lambdas;
           
            l_max = 1/max(lambdas);
            l_min = 1/min(lambdas);
            
            lambda = l_min + rand * (l_max - l_min);
            x2 = x + l_max * v;
            x = x + lambda * v;
            
        end
        p(:,i) = x2;
    end
end

