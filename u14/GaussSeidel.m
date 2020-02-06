function [u, uk] = GaussSeidel(A, b, u0, tol, uexact)
  U = triu(A, 1);
  Lstar = tril(A);
  Linv = inv(Lstar);
  uk = zeros(length(b), 1);
  uk(:, 1) = u0;
  k = 0;
  u = u0;
  while sqrt(dot(A*(u-uexact), (u-uexact))) >= tol
    u = Linv*(b - U*u);
    k += 1;
    uk(:, k+1) = u;
  endwhile
endfunction
