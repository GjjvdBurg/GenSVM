clear;
more off;

rand('state', 219038);

n = 6;
m = 5;

tmp = rand(n);
A = tmp + tmp' + n*eye(n);

for i=1:size(A, 1)
  for j=1:size(A, 2)
    if j >= i % only print the upper part
      fprintf('matrix_set(A, n, %i, %i, %.16f);\n', i-1, j-1, A(i, j));
    end
  end
end
fprintf('\n\n');

B = rand(n, m);
%b = vec(B); % this gives B in column-major order
%for i=1:numel(b)
%  fprintf('B[%i] = %.16f;\n', i-1, b(i));
%end

%for i=1:size(B, 1)
%  for j=1:size(B, 2)
%    fprintf('matrix_set(B, m, %i, %i, %.16f);\n', i-1, j-1, B(i, j));
%  end
%end
%fprintf('\n\n');
X = A \ B;

x = vec(X);
for i=1:numel(x)
  fprintf('mu_assert(fabs(B[%i] - %.16f) < 1e-14,\n"Incorrect value of B at %i");\n', i-1, x(i), i-1);
end
%
%for i=1:size(X, 1)
%  for j=1:size(X, 2)
%    fprintf('mu_assert(fabs(matrix_get(B, m, %i, %i) -\n%.16f) < 1e-14,\n"Incorrect value of B at %i, %i");\n', i-1, j-1, X(i, j), i-1, j-1);
%  end
%end
fprintf('\n\n');



%
%dposv(
%  UPLO, 'L'
%  N,    'n'
%  NRHS, 'm'
%  A,    'A'
%  LDA,  'n'
%  B,    'B'
%  LDB,  'n'
%  INFO)