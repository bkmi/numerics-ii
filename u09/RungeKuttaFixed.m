function [x, t] = RungeKuttaFixed(f, x0, I, tau, b, A, TOL)

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
  F = @(xx, zz)reshape(transpose(A*transpose(f(xx+zz))), s*d, 1);

  for i=2:N
    z = zeros(d, s);
    z_last = ones(d,s) .* (TOL + 1)^2;
    
    while norm(z - z_last) >= TOL
      z_last = z;
      z = tau * F(x(:,i-1), z);
    endwhile
    
    k(:,:,i) = f(x(:,i-1) + z);
    x(:,i) = x(:,i-1) + tau*k(:,:,i)*b;
  endfor  
endfunction