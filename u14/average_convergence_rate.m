function [rho] = average_convergence_rate(mat, uk, uexact)
error = @(M, x, xexact)sqrt(dot(M*(x-xexact),x-xexact));

u0 = uk(:, 1);
denom = error(mat, u0, uexact);
rho = 0;
for i = 2:size(uk, 2)
    numer = error(mat, uk(:, i), uexact);
    rho += (numer / denom)^(1/(i-1));
endfor
rho /= size(uk, 2);
end

