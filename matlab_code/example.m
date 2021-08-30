d = 2;
m = 7;

N = 2000;

A = randn(m, d);
b = rand(m, 1);

[x0, r]=get_cheb(A, b);

[A2, b2, center, radius, T, T_shift, t_r] = round_max_ellipsoid(A, b, x0);

H = inv(T);
%T2 = H'*H;

%H2 = inv(T2);
L1 = 2* sqrt(max(eig(H)));
L2 = 2*d*r;

X = randEllipsoid(d, N, H) + repmat(T_shift, [1 N]);

BX = boundary_hnr(A, b ,x0, 10000, 5);

[samples, avg_rho, acceptance_prob] = BilliardWalk_ellipsoid(A, b, x0, N, 17, L1, H);
[samples2, avg_rho2] = BilliardWalk(A, b, x0, N, 1, L2);



figure
plot(X(1,:), X(2,:), 'k.')
hold on
%figure
plot(BX(1,:), BX(2,:), 'k.')
%hold on
plot(samples(1,:), samples(2,:), 'r.');

