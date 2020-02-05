function A = create_laplacian (n)
  h = 1 / (1 + n);
  An = @(n) diag(-1 * ones(n-1,1), 1) ...
            + diag(-1 * ones(n-1,1), -1) ...
            + diag(4 * ones(n, 1));
            
  ACell = repmat({An(n)}, 1, n);
  A = blkdiag(ACell{:});
  A += diag(repmat(1, [1, (n-1)*n]), n);
  A += diag(repmat(1, [1, (n-1)*n]), -n);
endfunction