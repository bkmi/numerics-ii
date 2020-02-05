f = @(x, y) (x - x^2) * (y - y^2);
h_fn = @(k) 1/(k+1);
A_fn = @(k) h_fn(k)^(-2) * create_laplacian(k);
b_fn = @(k, grid) h_fn(k)^2 * create_rhs(k, f, grid);
error = @(M, x, xexact) sqrt(dot(M*(x-xexact),x-xexact));

tol = 10^-8;

n = 3;
A = A_fn(n);
grid = linspace(0, 1, n);
b = b_fn(n, grid);
u0 = zeros(n^2, 1);
uexact = A \ b;

[u, uk] = cg(A, b, u0, tol, uexact);

% [u_pre, uk_pre] = pcg_ben(A, b, u0, tol, uexact, pre);

n = 10;
A = A_fn(n);
grid = linspace(0, 1, n);
b = b_fn(n, grid);
u0 = zeros(n^2, 1);
uexact = A \ b;

[u, uk] = cg(A, b, u0, tol, uexact);

steps = size(uk, 2);
errors = zeros(steps, 1);
for i = 1:steps
  errors(i) = error(A, uk(:, i), uexact);
endfor

figure;
hold on;
plot(0:steps-1, errors);
title(strcat('CG n=', num2str(n)))
xlabel('steps');
ylabel('error');

% [u_pre, uk_pre] = pcg_ben(A, b, u0, tol, uexact, pre);

