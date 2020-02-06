f = @(x, y) (x - x^2) * (y - y^2);
h_fn = @(k) 1/(k+1);
A_fn = @(k) h_fn(k)^(-2) * create_laplacian(k);
U_fn = @(k, grid) h_fn(k)^2 * create_U(k, f, grid);
error = @(M, x, xexact)sqrt(dot(M*(x-xexact),x-xexact));
sgs = @(mat) tril(mat) * inv(diag(diag(mat))) * triu(mat);

cond = @(x)eig(x)(end)/eig(x)(1);

tol = 10^-8;

% cg
figure;
hold on;
labels = cell;

for n = 3:15;
  A = A_fn(n);
  grid = linspace(0, 1, n);
  uexact = U_fn(n, grid);
  b = A*uexact;
  u0 = zeros(n^2, 1);
  
  [u, uk] = cg(A, b, u0, tol, uexact);
  
  steps = size(uk, 2);
  errors = zeros(steps, 1);
  for i = 1:steps
    errors(i) = error(A, uk(:, i), uexact);
  endfor
  
  rho = average_convergence_rate(A, uk, uexact);
  
  plot(0:steps-1, errors);
  labels{end+1,1} = strcat('n=', num2str(n), ', rate = ', num2str(rho));
endfor

legend(labels)
title('Error vs Step for Conjugate Gradient Method')
xlabel('steps');
ylabel('error');

% pre cg
figure;
hold on;
labels = cell;
for n = 3:15;
  A = A_fn(n);
  grid = linspace(0, 1, n);
  uexact = U_fn(n, grid);
  b = A*uexact;
  u0 = zeros(n^2, 1);
  
  Binv = inv(sgs(A));
  pre = @(x) Binv * x;
  [u, uk] = pcg_ben(A, b, u0, tol, uexact, pre);
  % [u, uk] = cg(A, b, u0, tol, uexact);
  
  steps = size(uk, 2);
  errors = zeros(steps, 1);
  for i = 1:steps
    errors(i) = error(A, uk(:, i), uexact);
  endfor
  
  rho = average_convergence_rate(A, uk, uexact);
  
  plot(0:steps-1, errors);
  labels{end+1,1} = strcat('n=', num2str(n), ', rate = ', num2str(rho));
endfor

legend(labels)
title('Error vs Step for Preconditioned Conjugate Gradient Method')
xlabel('steps');
ylabel('error');