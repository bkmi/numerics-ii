c = [1; 2; 3];

## x is a column
f = @(x)[
  -c(1).*x(1,:)+c(3).*x(2,:).*x(3,:), 
  c(1).*x(1,:) - c(2).*x(2,:).^2 - c(3).*x(2,:).*x(3,:), 
  c(2).*x(2,:).^2
];
  
Df = @(x)[
  -c(1), c(3)*x(3), c(3)*x(2);
  c(1), -2*c(2)*x(2) - c(3)*x(3), -c(3)*x(2);
  0, 2*c(2)*x(2), 0
];


I = [0, 10];
tau_expos = 0:5;
taus = 0.01 * 2.^tau_expos

## For x0 = 0, 0, 0
# Implicit Euler
e_newton = zeros(size(taus));
e_fixed = zeros(size(taus));
e_simple = zeros(size(taus));

x0 = [0;0;0];
A = [[1]];
b = [1];
tol = tau^2;

for ind = 1:length(taus)
  tau = taus(ind)
  [x_newton, t] = RungeKuttaNewton(f, Df, x0, I, tau, b, A, tol);
  [x_fixed, t] = RungeKuttaFixed(f, x0, I, tau, b, A, tol);
  [x_simple, t] = RungeKuttaNewtonSimple(f, Df, x0, I, tau, b, A, tol);
  [_, x_truth] =  ode15s(@(t,x)f(x), t, x0);
  x_truth = transpose(x_truth);

  e_newton(ind) = mean(vecnorm(x_newton - x_truth));
  e_fixed(ind) = mean(vecnorm(x_fixed - x_truth));
  e_simple(ind) = mean(vecnorm(x_simple - x_truth));
endfor

figure
hold on
plot(taus, e_newton)
plot(taus, e_fixed)
plot(taus, e_simple)
legend('newton', 'fixed', 'simple')
xlabel('tau')
ylabel('mean square error')
title('Implicit Euler')
hold off

# Implicit Midpoint Rule
e_newton = zeros(size(taus));
e_fixed = zeros(size(taus));
e_simple = zeros(size(taus));

A = [[0.5]];
b = [1];
tol = tau^2;

for ind = 1:length(taus)
  tau = taus(ind)
  [x_newton, t] = RungeKuttaNewton(f, Df, x0, I, tau, b, A, tol);
  [x_fixed, t] = RungeKuttaFixed(f, x0, I, tau, b, A, tol);
  [x_simple, t] = RungeKuttaNewtonSimple(f, Df, x0, I, tau, b, A, tol);
  [_, x_truth] =  ode15s(@(t,x)f(x), t, x0);
  x_truth = transpose(x_truth);

  e_newton(ind) = mean(vecnorm(x_newton - x_truth));
  e_fixed(ind) = mean(vecnorm(x_fixed - x_truth));
  e_simple(ind) = mean(vecnorm(x_simple - x_truth));
endfor

figure
hold on
plot(taus, e_newton)
plot(taus, e_fixed)
plot(taus, e_simple)
legend('newton', 'fixed', 'simple')
xlabel('tau')
ylabel('mean square error')
title('Implicit Midpoint')
hold off

## For x0 = 1, 0, 0
# Implicit Euler
e_newton = zeros(size(taus));
e_fixed = zeros(size(taus));
e_simple = zeros(size(taus));

x0 = [1;0;0];
A = [[1]];
b = [1];
tol = tau^2;

for ind = 1:length(taus)
  tau = taus(ind)
  [x_newton, t] = RungeKuttaNewton(f, Df, x0, I, tau, b, A, tol);
  [x_fixed, t] = RungeKuttaFixed(f, x0, I, tau, b, A, tol);
  [x_simple, t] = RungeKuttaNewtonSimple(f, Df, x0, I, tau, b, A, tol);
  [_, x_truth] =  ode15s(@(t,x)f(x), t, x0);
  x_truth = transpose(x_truth);

  e_newton(ind) = mean(vecnorm(x_newton - x_truth));
  e_fixed(ind) = mean(vecnorm(x_fixed - x_truth));
  e_simple(ind) = mean(vecnorm(x_simple - x_truth));
endfor

figure
hold on
plot(taus, e_newton)
plot(taus, e_fixed)
plot(taus, e_simple)
legend('newton', 'fixed', 'simple')
xlabel('tau')
ylabel('mean square error')
title('Implicit Euler')
hold off

# Implicit Midpoint Rule
e_newton = zeros(size(taus));
e_fixed = zeros(size(taus));
e_simple = zeros(size(taus));

A = [[0.5]];
b = [1];
tol = tau^2;

for ind = 1:length(taus)
  tau = taus(ind)
  [x_newton, t] = RungeKuttaNewton(f, Df, x0, I, tau, b, A, tol);
  [x_fixed, t] = RungeKuttaFixed(f, x0, I, tau, b, A, tol);
  [x_simple, t] = RungeKuttaNewtonSimple(f, Df, x0, I, tau, b, A, tol);
  [_, x_truth] =  ode15s(@(t,x)f(x), t, x0);
  x_truth = transpose(x_truth);

  e_newton(ind) = mean(vecnorm(x_newton - x_truth));
  e_fixed(ind) = mean(vecnorm(x_fixed - x_truth));
  e_simple(ind) = mean(vecnorm(x_simple - x_truth));
endfor

figure
hold on
plot(taus, e_newton)
plot(taus, e_fixed)
plot(taus, e_simple)
legend('newton', 'fixed', 'simple')
xlabel('tau')
ylabel('mean square error')
title('Implicit Midpoint')
hold off