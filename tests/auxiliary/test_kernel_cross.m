
function test_kernel_cross(kerneltype)

clc;
rand('state', 123456);
more off;

n_1 = 10;
n_2 = 5;
m = 3;

X_1 = rand(n_1, m);
Z_1 = [ones(n_1, 1), X_1];
X_2 = rand(n_2, m);
Z_2 = [ones(n_2, 1), X_2];

for ii=1:n_1
  for jj=1:m+1
    fprintf("matrix_set(data_1->RAW, data_1->m+1, %i, %i, %.16f);\n", ii-1, jj-1, Z_1(ii, jj));
  end
end

fprintf('\n');
for ii=1:n_2
  for jj=1:m+1
    fprintf("matrix_set(data_2->RAW, data_2->m+1, %i, %i, %.16f);\n", ii-1, jj-1, Z_2(ii, jj));
  end
end


K = zeros(n_2, n_1);
if strcmp(kerneltype, 'poly')
  # Polynomial kernel
  # (gamma * <x_1, x_2> + c)^d
  gamma = 1.5;
  c = 3.0;
  d = 1.78;
  
  for ii=1:n_2
    for jj=1:n_1
      K(ii, jj) = (gamma * (X_2(ii, :) * X_1(jj, :)') + c)^d;
    end
  end
elseif strcmp(kerneltype, 'rbf')
  # RBF kernel
  # exp(-gamma * norm(x1 - x2)^2)
  gamma = 0.348
  for ii=1:n_2
    for jj=1:n_1
      K(ii, jj) = exp(-gamma * sum((X_2(ii, :) - X_1(jj, :)).^2));
    end
  end
elseif strcmp(kerneltype, 'sigmoid')
  # Sigmoid kernel
  # tanh(gamma * <x_1, x_2> + c)
  gamma = 1.23;
  c = 1.6;
  for ii=1:n_2
    for jj=1:n_1
      K(ii, jj) = tanh(gamma * (X_2(ii, :) * X_1(jj, :)') + c);
    end
  end
end

fprintf('\n');
for ii=1:n_2
  for jj=1:n_1
    fprintf("mu_assert(fabs(matrix_get(K2, data_1->n, %i, %i) -\n %.16f) < eps,\n\"Incorrect K2 at %i, %i\");\n", ii-1, jj-1, K(ii, jj), ii-1, jj-1);
  end
end