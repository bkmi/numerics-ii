function b = create_rhs (k, f, grid)
  assert(length(grid) == k)
  b = zeros(k^2, 1);
  for i = 1:k
    for j = 1:k
      b(k * (i-1) + j) = f(grid(i), grid(j));
    endfor
  endfor

endfunction
