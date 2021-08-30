function [E, acceptance_prob] = determine_covariance(A, b, x0, N, W, L, sigma, prob_max, prob_min)

    [Q, D] = eig(sigma);
    n = size(sigma, 2);
    
    max_eig = max(diag(D));
    min_eig = min(diag(D));
    
    identity = eye(n);
    
    e = 1.001;
    max_scale = (max_eig - min_eig * e) / (e - 1); % this factor corresponds to the classical billiard walk
    min_scale = 0;
    
    while(true)
        
        med_scale = (max_scale + min_scale) / 2;
        Dnew = D + med_scale * identity;
        
        max_eig = max(diag(Dnew))
        min_eig = min(diag(Dnew))
        
        E = Q * Dnew * Q';
        [~, ~, acceptance_prob] = BilliardWalk_ellipsoid(A, b, x0, N, W, L, E);
        acceptance_prob
        
        if (acceptance_prob < prob_max && acceptance_prob > prob_min)
            return
        elseif (acceptance_prob > prob_max)
            max_scale = med_scale;
        else
            min_scale = med_scale;
        end
    end

end
