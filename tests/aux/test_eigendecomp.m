clear;
rand('state', 123456);
more off; # for Octave

n = 10;

A = rand(n, n);
K = A'*A;


for ii=1:n
  for jj=1:n
    fprintf("matrix_set(K, n, %i, %i, %.16f);\n", ii-1, jj-1, K(ii, jj));
  end
end

[P, Sigma] = eig(K);


eigenvalues = diag(Sigma);
ratios = eigenvalues ./ eigenvalues(end, end);

cutoff = 1e-2;

realP = fliplr(P(:, ratios > cutoff));
realSigma = sqrt(flipud(eigenvalues(ratios > cutoff)));

r = sum(ratios > cutoff);

assert(r == size(realP, 2));

fprintf('\n');
for ii=1:n
  for jj=1:r
    fprintf("mu_assert(fabs(fabs(matrix_get(P, r, %i, %i)) -\nfabs(%.16f)) < eps,\n\"Incorrect P at %i, %i\");\n", ii-1, jj-1, realP(ii, jj), ii-1, jj-1);
  end
end

fprintf('\n');
for jj=1:r
  fprintf("mu_assert(fabs(Sigma[%i] - %.16f) < eps,\n\"Incorrect Sigma at %i\");\n", jj-1, realSigma(jj), jj-1);
end

