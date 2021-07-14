A = polytope.A;
b = polytope.b;

d = size(A, 2);

normA = sqrt(sum(A.^2, 2));
A = A ./ repmat(normA, [1, d]);
b = b ./ normA;

x = polytope.center;
radius = polytope.radius;
L = 30 * d * radius;

tic
[samples1, avg_rho1, acceptance_prob] = BilliardWalk_hessian(A, b, x, 30000, 1, L);
tim1 = toc;

tic
[samples2, avg_rho2] = BilliardWalk(A, b, x, 30000, 1, L);
tim2 = toc;

% compute the PSRF of each marginal
R1 = psrf(samples1');
R2 = psrf(samples2');
    
% compute the Effective Sample Size of each marginal
ess1 = effective_sample_size(samples1)';
ess2 = effective_sample_size(samples2)';

%plot(samples(1,:), samples(2,:), '.')