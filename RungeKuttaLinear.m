function  [x, t, k] = RungeKuttaLinear (M, x0, I, tau, b, A)
  s = length(b);
  d = length(x0);
  n = (I(2) - I(1))/tau;
  c = sum(A, 2);
  t = I(1) + tau * c;
  t = t'
  
  x = zeros(d, n)
  t = 
  for i = 1:n
    
  soln = fsolve(system(M, x0, tau, b, A), zeros(6,1));
  x = soln(1:2);
  k = soln(3:6);
endfunction

function f = system (M, xn, tau, b, A)
  y = @(xnext, k1, k2) xn + tau * ( b(1) * k1 + b(2) * k2) - xnext;
  kone = @(k1, k2) M * (xn + tau * (A(1,1) * k1 + A(1,2) * k2)) - k1;
  ktwo = @(k1, k2) M * (xn + tau * (A(2,1) * k1 + A(2,2) * k2)) - k2;
  f = @(x)[y(x(1:2), x(3:4), x(5:6)),
           kone(x(3:4), x(5:6)),
           ktwo(x(3:4), x(5:6))];
endfunction