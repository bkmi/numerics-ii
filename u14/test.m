a = [1, -0.1; 0.25, 1.1];
A = a' * a;
b = [1; 2];
u0 = [0.5; 0.4];
tol = 0.01;

uexact = A \ b;

[u, uk] = cg(A, b, u0, tol, uexact);

% B = [1, 0; 0, 1];
B = [1, 0.1; 0.1, 1];
Binv = inv(B);
pre = @(x)Binv * x;
[u_pre, uk_pre] = pcg_ben(A, b, u0, tol, uexact, pre);