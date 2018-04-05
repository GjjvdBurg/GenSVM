function test_kernel_pre()

  kerneltype = 'rbf';
  rand('state', 123456);
  n = 10;
  m = 5;
  cutoff = 5e-3;

  X = rand(n, m);
  Z = [ones(n, 1), X];

  set_matrix(Z, "data->Z", "data->m+1");

  K = zeros(n, n);
  if strcmp(kerneltype, 'poly')
    % Polynomial kernel
    % (gamma * <x_1, x_2> + c)^d
    gamma = 1.5;
    c = 3.0;
    d = 1.78;

    for ii=1:n
      for jj=1:n
        K(ii, jj) = (gamma * (X(ii, :) * X(jj, :)') + c)^d;
      end
    end
  elseif strcmp(kerneltype, 'rbf')
    % RBF kernel
    % exp(-gamma * norm(x1 - x2)^2)
    gamma = 0.348
    for ii=1:n
      for jj=1:n
        K(ii, jj) = exp(-gamma * sum((X(ii, :) - X(jj, :)).^2));
      end
    end
  elseif strcmp(kerneltype, 'sigmoid')
    % Sigmoid kernel
    % tanh(gamma * <x_1, x_2> + c)
    gamma = 1.23;
    c = 1.6;
    for ii=1:n
      for jj=1:n
        K(ii, jj) = tanh(gamma * (X(ii, :) * X(jj, :)') + c);
      end
    end
  end

  [P, values] = eig(K);

  eigenvalues = diag(values);
  ratios = eigenvalues ./ eigenvalues(end, end);

  realP = fliplr(P(:, ratios > cutoff));
  realSigma = sqrt(flipud(eigenvalues(ratios > cutoff)));

  assert_matrix(realSigma, "data->Sigma", "1");

  r = sum(ratios > cutoff);
  fprintf("mu_assert(data->r == %i);\n", r);

  M = realP * diag(realSigma);

  newZ = [ones(n, 1) M];
  assert_matrix_abs(newZ, "data->Z", "data->r+1");

  assert_matrix(Z, "data->RAW", "data->m+1");

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
