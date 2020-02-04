a = [1, -0.1; 0.25, 1.1];
A = a' * a;
b = [1; 2];
u0 = [0.5; 0.4];
tol = 0.01;

uexact = A \ b;

[u, uk] = cg(A, b, u0, tol, uexact);