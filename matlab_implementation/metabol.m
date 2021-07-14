S = [1 2 -3];
n=3;
m=2;
A = [eye(n); -eye(n)];
b = [2; 3; 5; 1.5; 1; 3];
N = null(S);
AA = A*N;
x = [0;0];
AA*x-b
[p, avg_rho] = BilliardWalk(AA, b, x, 1000, 3, 3);

X = N*p;
plot3(X(1,:), X(2,:), X(3,:),'.')

grid on

