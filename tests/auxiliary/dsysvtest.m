more off;
clear;

rand('state', 891716);

n = 6;
m = 5;


tmp = rand(n);
% A is symmetric, but not necessarily p.s.d.
A = tmp + tmp';
clear tmp;

for i=1:size(A, 1)
  for j=1:size(A, 2)
    if j >= i % only print the upper part
      fprintf('matrix_set(A, n, %i, %i, %.16f);\n', i-1, j-1, A(i, j));
    end
  end
end
fprintf('\n\n');

B = rand(n, m);
b = vec(B); % this gives B in column-major order
for i=1:numel(b)
  fprintf('B[%i] = %.16f;\n', i-1, b(i));
end

X = A \ B;

x = vec(X);
for i=1:numel(x)
  fprintf('mu_assert(fabs(B[%i] - %.16f) < 1e-14,\n"Incorrect value of B at %i");\n', i-1, x(i), i-1);
end

fprintf('\n\n');

