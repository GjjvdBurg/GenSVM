
function kernel_tests(kerneltype)

rand('state', 123456);
more off;

n = 10;
m = 3;

X = rand(n, m);
Z = [ones(n, 1), X];

for ii=1:n
  for jj=1:m+1
    fprintf("matrix_set(data->RAW, data->m+1, %i, %i, %.16f);\n", ii-1, jj-1, Z(ii, jj));
  end
end


K = zeros(n, n);
if strcmp(kerneltype, 'poly')
  # Polynomial kernel
  # (gamma * <x_1, x_2> + c)^d
  gamma = 1.5;
  c = 3.0;
  d = 1.78;
  
  for ii=1:n
    for jj=1:n
      K(ii, jj) = (gamma * (X(ii, :) * X(jj, :)') + c)^d;
    end
  end
elseif strcmp(kerneltype, 'rbf')
  # RBF kernel
  # exp(-gamma * norm(x1 - x2)^2)
  gamma = 0.348
  for ii=1:n
    for jj=1:n
      K(ii, jj) = exp(-gamma * sum((X(ii, :) - X(jj, :)).^2));
    end
  end
elseif strcmp(kerneltype, 'sigmoid')
  # Sigmoid kernel
  # tanh(gamma * <x_1, x_2> + c)
  gamma = 1.23;
  c = 1.6;
  for ii=1:n
    for jj=1:n
      K(ii, jj) = tanh(gamma * (X(ii, :) * X(jj, :)') + c);
    end
  end
end

fprintf('\n');
for ii=1:n
  for jj=1:n
    fprintf("mu_assert(fabs(matrix_get(K, data->n, %i, %i) -\n %.16f) < eps,\n\"Incorrect K at %i, %i\");\n", ii-1, jj-1, K(ii, jj), ii-1, jj-1);
  end
end