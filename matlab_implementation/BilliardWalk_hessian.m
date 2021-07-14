function [samples, avg_rho, acceptance_prob] = BilliardWalk_hessian(A, b, x0, N, W, L)
% Billiard Walk to sample from the uniform distribution
        
    d = size(A, 2);
    samples = zeros(d, N);
    mu = zeros(d, 1);
    
    h1 = 0;
    h2 = 0;
    
    avg_rho = 0;
    acceptance_prob = 0;
    
    h = waitbar(0,'Computing samples...');
    for i = 1:N
        for j = 1:W
            x = x0;
            T = -log(rand) * L;
            D = diag(b - A* x);
            Hessian = A' * D^2 * A;
            v = mvnrnd(mu, inv(Hessian), 1)';
            T = T/norm(v);
            h1 = 0.5 * log(Hessian) + 0.5 * (v' * Hessian * v);
            
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
            D = diag(b - A* x);
            Hessian = A' * D^2 * A;
            h2 = 0.5 * log(Hessian) + 0.5 * (v' * Hessian * v);
            
            log_prob = log(rand);
            if (log_prob < (h2 - h1))
                x0 = x;
                acceptance_prob = acceptance_prob + 1;
            end
            avg_rho = avg_rho + rho;
        end
        
        samples(:, i) = x0;
        waitbar(i/N);
    end
    avg_rho = avg_rho / (N * W);
    acceptance_prob = acceptance_prob / (N * W);
    close(h);
end

