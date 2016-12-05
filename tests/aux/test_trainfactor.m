clc;
more off; # Octave
rand('state', 123456);

n = 10;
m = 5;
K = 3;
r = 7;

P = rand(n, r);
Sigma = rand(r, 1);

for ii=1:n
  for jj=1:r
    fprintf("matrix_set(P, r, %i, %i, %.16f);\n", ii-1, jj-1, P(ii, jj));
  end
end

fprintf('\n');
for ii=1:r
  fprintf("Sigma[%i] = %.16f;\n", ii-1, Sigma(ii));
end

Z = [ones(n, 1), P*diag(Sigma)];

for ii=1:n
  for jj=1:r+1
    fprintf("mu_assert(fabs(matrix_get(Z, r+1, %i, %i) -\n%.16f) < eps,\n\"Incorrect Z at %i, %i\");\n", ii-1, jj-1, Z(ii, jj), ii-1, jj-1);
  end
end



  