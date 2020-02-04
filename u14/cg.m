function [u, uk] = cg(A, b, u0, tol, uexact)
    n = length(b);
    uk = zeros(n, n); % At least the minimum size for cg.
    uk(:, 1) = u0;
    
    r0 = b - A * u0;
    e0 = r0;
    k = 0;
    
    u = u0;
    r = r0;
    e = e0;
    r_norm_old = r' * r;
    while norm(u - uexact) >= tol
      Ae = A * e;
      alpha = r_norm_old / (Ae' * e);
      
      u = u + alpha * e;
      r = r - alpha * Ae;
      r_norm_new = r' * r;
      beta = r_norm_new / r_norm_old;
      e = r + beta * e;
      r_norm_old = r_norm_new
      
      uk(:, k + 1) = u;
      k += 1;
    endwhile
endfunction