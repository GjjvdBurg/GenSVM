function test_testfactor()
  rand('state', 123456);
  more off; # Octave
  
  train_n = 10;
  train_r = 3;
  test_n = 5;
  
  
  P = rand(train_n, train_r);
  K2 = rand(test_n, train_n);
  Sigma = diag(rand(1, train_r));
  
  train_Z = [ones(train_n, 1), P];
  
  set_matrix(train_Z, "train->Z", "train->r+1");
  set_matrix(diag(Sigma), "Sigma", "1");
  set_matrix(K2, "K2", "train->n");
  
  test_factor = K2 * P * Sigma^(-2);
  test_Z = [ones(test_n, 1) test_factor];
  assert_matrix(test->Z, "test->Z", "test->r+1");
  
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
