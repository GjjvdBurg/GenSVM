function [V] = test_train_kernel()
  
  clear;
  more off;
  rand('state', 654321);
  
  n = 10;
  m = 5;
  classes = 4;
  cutoff = 5e-3;
  
  X = rand(n, m);
  Z = [ones(n, 1), X];
  set_matrix(Z, "data->Z", "data->m+1");
  
  y = [2 1 3 2 3 2 4 1 3 4];
  set_matrix(y, "data->y", "1");
  
  p = 1.2143;
  kappa = 0.90298;
  lambda = 0.00219038;
  epsilon = 1e-15;
  rho = ones(n, 1);

  K = zeros(n, n);
  # RBF kernel
  # exp(-gamma * norm(x1 - x2)^2)
  gamma = 0.348
  for ii=1:n
    for jj=1:n
      K(ii, jj) = exp(-gamma * sum((X(ii, :) - X(jj, :)).^2));
    end
  end

  [P, Sigma] = eig(K);

  eigenvalues = diag(Sigma);
  ratios = eigenvalues ./ eigenvalues(end, end);
 
  realP = fliplr(P(:, ratios > cutoff));
  realSigma = flipud(eigenvalues(ratios > cutoff));
  
  assert_matrix(realSigma, "data->Sigma", "1");
  
  r = sum(ratios > cutoff);
  fprintf("mu_assert(data->r == %i);\n", r);
  
  M = realP * diag(realSigma);
  size(M)
  
  assert_matrix(Z, "data->RAW", "data->m+1");
  
  seedV = zeros(size(M, 2) + 1, classes - 1);  
  [W, t] = msvmmaj(M, y, rho, p, kappa, lambda, epsilon, 'show', 0, seedV);
  V = [t'; W];
  
  fprintf('\n');
  assert_matrix_abs(V, "model->V", "model->K-1");

end

function set_matrix(A, name, cols)
  for ii=1:size(A, 1)
    for jj=1:size(A, 2)
      fprintf("matrix_set(%s, %s, %i, %i, %.16f);\n", name, cols, ii-1, jj-1, A(ii, jj));
    end
  end
  fprintf("\n");
end

function assert_matrix(A, name, cols)
  for ii=1:size(A, 1)
    for jj=1:size(A, 2)
          fprintf(["mu_assert(fabs(matrix_get(%s, %s, %i, %i) -\n%.16f) <", ...
                   " eps,\n\"Incorrect %s at %i, %i\");\n"], name, cols, ...
                   ii-1, jj-1, A(ii, jj), name, ii-1, jj-1);
    end
  end
  fprintf("\n");
end


function assert_matrix_abs(A, name, cols)
  for ii=1:size(A, 1)
    for jj=1:size(A, 2)
          fprintf(["mu_assert(fabs(fabs(matrix_get(%s, %s, %i, %i)) -\nfabs(%.16f)) <", ...
                   " eps,\n\"Incorrect %s at %i, %i\");\n"], name, cols, ...
                   ii-1, jj-1, A(ii, jj), name, ii-1, jj-1);
    end
  end
  fprintf("\n");
end