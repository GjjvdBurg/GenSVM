
function test_kernel_post()
  more off;
  rand('state', 123456);
  
  n_tr = 10;
  n_tt = 8;
  m = 5;
  K = 3;
  
  X_tr = rand(n_tr, m);
  Z_tr = [ones(n_tr, 1), X_tr];
  set_matrix(Z_tr, "train->RAW", "train->r+1");
  fprintf("train->Z = train->RAW;\n");
  
  X_tt = rand(n_tt, m);
  Z_tt = [ones(n_tt, 1), X_tt];
  set_matrix(Z_tt, "test->RAW", "test->m+1");
  fprintf("test->Z = test->RAW;\n");
  
  gamma = 1.132;
  K2 = zeros(n_tt, n_tr);
  for ii=1:n_tt
    for jj=1:n_tr
      K2(ii, jj) = exp(-gamma * sum((X_tt(ii, :) - X_tr(jj, :)).^2));
    end
  end
  
  M = X_tr;
  
  Sigma = rand(m, 1);
  set_matrix(Sigma, "train->Sigma", "1");
  
  N = K2 * M * diag(Sigma)^(-2);
  test_Z = [ones(size(N, 1), 1), N];
  
  assert_matrix(test_Z, "test->Z", "test->r+1");
  
  
  

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