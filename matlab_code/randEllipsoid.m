function X = randEllipsoid(d, N, sigma)
    
    H = inv(sigma);
    L = chol(H)';
    X = randsphere(N,d,1)';
    
    X = L * X;
end

