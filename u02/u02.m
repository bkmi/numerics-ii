twostep = @(xk, xkminus1, l, tau) xkminus1 - 2 * tau * l * xk;
solution = @(t, l, x0) x0 * exp(-l*t);

lamb = 1.0;
x0 = 1.0;
tmin = 0.0;
tmax = 1.0;

hold on

grid = linspace(tmin, tmax, 1000);
y_soln = solution(grid, lamb, x0);
plot(grid, y_soln)

for n = [10, 100, 1000];
    grid = linspace(tmin, tmax, n);
    tau = 1/n;
    x1 = x0 - tau * lamb * x0;
    y_twostep = [x0 x1];
    for tk = grid(3:end)
        y_twostep(end+1) = twostep(y_twostep(end), y_twostep(end-1), lamb, tau);
    end
    plot(grid, y_twostep)
end
