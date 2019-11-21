function  [x, t, k] = RungeKuttaLinear (M, x0, I, tau, b, A)
  s = length(b);
  d = length(x0);
  n = (I(2) - I(1))/tau;
  
  x = zeros(d, n);
  t = [I(1), zeros(1, n - 1)];
  k = zeros(n, d, s);
  
  for i = 2:n
    t(i) = I(1) + tau * i;
    soln = fsolve(system(M, x0, tau, b, A), [x0; ones(4,1)]);
    x(:, i) = soln(1:2);
    k(i, :, 1) = soln(3:4);
    k(i, :, 2) = soln(5:6);
  end
endfunction

function f = system (M, xn, tau, b, A)
  y = @(xnext, k1, k2) xn + tau * ( b(1) * k1 + b(2) * k2) - xnext;
  kone = @(k1, k2) M * (xn + tau * (A(1,1) * k1 + A(1,2) * k2)) - k1;
  ktwo = @(k1, k2) M * (xn + tau * (A(2,1) * k1 + A(2,2) * k2)) - k2;
  f = @(x)[y(x(1:2), x(3:4), x(5:6)),
           kone(x(3:4), x(5:6)),
           ktwo(x(3:4), x(5:6))];
endfunction