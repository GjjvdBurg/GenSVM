function [V] = test_optimize()

clear;
more off;
rand('state', 902183);

n = 8;
m = 3;
K = 4;

X = rand(n, m);
set_matrix(X);

y = [2 1 3 2 3 2 4 1];

set_vector(y);

p = 1.2143;
kappa = 0.90298;
lambda = 0.00219038;
epsilon = 1e-15;

rho = rand(n, 1);
set_vector(rho);

[W, t] = msvmmaj(X, y, rho, p, kappa, lambda, epsilon);

V = [t'; W];

assert_matrix(V, "model->V", "model->K-1");

end

function set_matrix(A)
  for i=1:size(A, 1)
    for j=1:size(A, 2)
      fprintf('matrix_set(A, %i, %i, %i, %.16f);\n', size(A, 2), i-1, j-1, A(i, j));
    end
  end
end

function set_vector(a)
  for i=1:numel(a)
    fprintf('A[%i] = %.16f;\n', i-1, a(i));
  end
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