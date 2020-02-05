h_fn = @(k) 1/(k+1);
A_fn = @(k) h_fn(k)^(-2) * create_laplacian(k);
f = @(x, y) (x - x^2) * (y - y^2);

n = 5;
A = A_fn(n);
grid = linspace(0, 1, n);

b = zeros(n^2, 1);
for i = 1:n
  for j = 1:n
    b(n * (i-1) + j) = f(grid(i), grid(j));
  endfor
endfor
b *= h_fn(n)^2;


A \ b