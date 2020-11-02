m=6; d=6;
X = zeros(2*m, 2*m);
Y = zeros(2*m, 2*m);

M = 2*rand(m) - 1;
matrices{1} = -M*M' - eye(m);

for i=1:d
    M = 2*rand(m/2) - 1;
    M = M + M';
    matrices{i+1} = blkdiag(M, -M);
end

x=zeros(d,0);
v= randn(d,1);
c = randn(d,1);
c=c/norm(c);
v=v/norm(v);

C = lmi_eval(matrices, x, true);
B = lmi_eval(matrices, v, false);
A = lmi_eval(matrices, c, false);

X(1:m,1:m) = B;
X(1:m, (m+1):(2*m)) = C;
X((m+1):(2*m), 1:m) = C;

Y(1:m,1:m) = -A;
Y((m+1):(2*m), (m+1):(2*m)) = C;

 [V,D] = eig(X,Y);
 
 diag(D)
 
R = zeros(2*m, 2*m);
P = zeros(2*m, 2*m);

R(1:m,1:m) = A;
R((m+1):(2*m), (m+1):(2*m)) = eye(m);

P(1:m,1:m) = -B;
P(1:m, (m+1):(2*m)) = -C;
P((m+1):(2*m), 1:m) = eye(m);


