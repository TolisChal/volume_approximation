d = 20;
m = 200;

N = 1000;

A = randn(m, d);
b = rand(m, 1);

[x0, r]=get_cheb(A, b);

[A2, b2, center, radius, T, T_shift, t_r] = round_max_ellipsoid(A, b, x0);

H = inv(T);
L = 2*d*r;

% try to achieve an acceptance probability in [0.65, 0.72]
[E, acceptance_probability] = determine_covariance(A, b, x0, N, 3, L, H, 0.72, 0.65);



