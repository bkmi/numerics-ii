function [p, q, t] = stoermer_verlet (m, g, r0, p0, q0, I, tau)

  d = length(p0);
  t = I(1):tau:I(2);
  N = length(t);
  dt = (I(2)-I(1))/N;
  p = zeros(d, N);
  p(:,1) = p0;
  q = zeros(d, N);
  q(:,1) = q0;
  
  f_dq = @(pp) pp/m;
  f_dp = @(qq) -m * g / r0 * cos(qq);
  step = @(q1, q2, dpp) 2 * q2 - q1 + dpp * dt^2;
  velo = @(q1, q2) (q2 - q1) / (2 * dt);
  
  for i=2:N
    dq = f_dq(p(:, i-1));
    dp = f_dp(q(:, i-1));
    if i == 2
      q(:, i) = q(:, i-1) + dq * dt + dp * dt^2 * 0.5;
    else
      q(:, i) = step(q(:, i-2), q(:, i-1), p(:, i-1));
    endif
    p(:, i) = velo(q(:, i-1), q(:, i)); 
  endfor
  
  
  

endfunction
