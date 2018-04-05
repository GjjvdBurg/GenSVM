function [V] = test_train()

clear;
more off;
rand('state', 123456);

n = 10;
m = 3;
K = 4;

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
seedV = rand(m+1, K-1);
set_matrix(seedV, "seed->V", "data->K-1");

[W, t] = msvmmaj(X, y, rho, p, kappa, lambda, epsilon, 'show', 0, seedV);
V = [t'; W];

fprintf('\n');
assert_matrix(V, "model->V", "model->K-1");

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