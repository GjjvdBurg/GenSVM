more off;
clear;

rand('state', 21353242);

n = 5;
m = 3;
K = 3;

p = 1.213;

y = [   2   1   3   2   3]';

for i=1:numel(y)
  fprintf('data->y[%i] = %.16f;\n', i-1, y(i));
end

H = 2*rand(n, m);

for i=1:size(H, 1)
  for j=1:size(H, 2)
    fprintf('matrix_set(A, %i, %i, %i, %.16f);\n', size(H, 2), i-1, j-1, H(i, j));
  end
end

R = zeros(n, K);
I = eye(K);
for i=1:n
  R(i, :) = I(y(i, :), :);
end
R = ~logical(R);

omega = (1/p)*(sum((H.^p).*R,2)).^(1/p - 1);

for i=1:n
  fprintf('mu_assert(fabs(gensvm_calculate_omega(model, %i) -\n%.16f) < 1e-14,\n"Incorrect omega at %i");\n', i-1, omega(i), i-1);
end