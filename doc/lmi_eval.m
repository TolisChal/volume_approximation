function [A] = lmi_eval(matrices, p, a0)
    

    d = length(p);
    m = size(matrices{1},1);
    A = zeros(m,m);
    for i=2:(d+1)
        A = A + matrices{i}*p(i-1);
    end
    
    if(a0)
        A = A + matrices{1};
    end
end