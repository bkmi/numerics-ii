function [x, t] = RungeKuttaNewton(f, Df, x0, I, tau, b, A, TOL)

  % reformulation:
  % k_i = f(g_i) ,  i=1,...s
  % z_i = g_i - x , i=1,...,s
  % Z = [z_1;...;z_s].

  b = b(:);
  d = length(x0);
  s = length(b);
  t = I(1):tau:I(2);
  N = length(t);
  x = zeros(d, N);
  k = zeros(d,s,N);
  x(:,1) = x0;
  DF = @(xx, zz)DF_helper(A, Df, xx, zz);
  F = @(xx, zz)reshape(transpose(A*transpose(f(xx+zz))), s*d, 1);

  for i=2:N  
    Z = zeros(s*d, 1);
    Z_last = ones(s*d, 1) .* (TOL + 1)^2;
    while norm(Z - Z_last) >= TOL
      z = reshape(Z, d, s);
      DG0 = eye(s*d) - tau*DF(x(:,i-1), z);
      Z = reshape(z, s*d, 1);
      G0 = -(Z - tau*F(x(:,i-1), z));
      DZ = DG0\G0;
      Z_last = Z;
      Z = Z + DZ;
    endwhile
    
    X = kron(ones(s,1), x(:,i-1));
    k(:,:,i) = f(reshape(X+Z, d, s));
    x(:,i) = x(:,i-1) + tau*k(:,:,i)*b;
  endfor  
endfunction

function out = DF_helper(A, Df, x, z)
  d = size(z)(1);
  s = size(z)(2);

  out = zeros(s*d, s*s);
  for i = 1:s
    for j = 1:s
      out((i-1)*d+1:(i-1)*d+d, (j-1)*d+1:(j-1)*d+d) = A(i,j).*Df(x + z(:, j));
    endfor
  endfor
endfunction
