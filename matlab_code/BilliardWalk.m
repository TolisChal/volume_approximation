function [samples, avg_rho] = BilliardWalk(A, b, x, N, W, L)
% Billiard Walk to sample from the uniform distribution
        
    d = size(A, 2);
    samples = zeros(d, N);
    avg_rho = 0;
    
    row_norms = sqrt(sum(A.^2,2));
    A = diag(1./row_norms)*A;
    b = diag(1./row_norms)*b;
    
    h = waitbar(0,'Computing samples...');
    for i = 1:N
        for j = 1:W
            
            T = -log(rand) * L;
            v = get_direction(d);
            
            rho = 0;
            while (true)
            
                % compute the intersection of the line x + l*v with the
                % boundary of the hypercube [-1, 1]^n
                lambdas =  A*v ./ (b - A*x);
                [l_max, pos_max] = max(lambdas);
                l_max = 1 / l_max;

                % pick a uniformly distributed point from the segment
                lambda = 0.995 * l_max;
                if (T <= lambda)
                    x = x + T * v;
                    break;
                end
                s = A(pos_max, :)';
                % update the current point of the random walk
                x = x + lambda * v;
                %reflect the ray
                v = v - (2*(v'*s))*s;
                rho = rho+1;
                T = T - lambda;
            end
            avg_rho = avg_rho + rho;
        end
        
        samples(:, i) = x;
        waitbar(i/N);
    end
    avg_rho = avg_rho / (N * W);
    close(h);
end