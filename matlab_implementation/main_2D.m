A = randn(7, 2);
b = rand(7, 1);

d = size(A, 2);

normA = sqrt(sum(A.^2, 2));
A = A ./ repmat(normA, [1, d]);
b = b ./ normA;

x = zeros(2, 1);
L = 4;

[samples, avg_rho] = BilliardWalk(A, b, x, 5000, 1, L);
%[samples, avg_rho] = BilliardWalk_hessian(A, b, x, 5000, 1, L);

% compute the PSRF of each marginal
R = psrf(samples');
    
% compute the Effective Sample Size of each marginal
ess = effective_sample_size(samples)';

plot(samples(1,:), samples(2,:), '.')