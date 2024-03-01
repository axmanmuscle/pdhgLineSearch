n = 5;

A = randn(n);
A = A / norm(A);
theta = 1/norm(A)^2 - 1e-6;

bt = chol((1/theta)*eye(n) - A*A');
B = bt';

A*A' + B*B'

norm(A)
norm(B)