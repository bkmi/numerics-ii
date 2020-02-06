function [u, uk] = Jacobi(A, b, u0, tol, uexact)
  D = diag(diag(A));
  Dinv = inv(D);
  R = A - D;
  uk = zeros(length(b), 1);
  uk(:, 1) = u0;
  k = 0;
  u = u0;
  while sqrt(dot(A*(u-uexact), (u-uexact))) >= tol
    u = Dinv*(b - R*u);
    k += 1;
    uk(:, k+1) = u;
  endwhile
endfunction